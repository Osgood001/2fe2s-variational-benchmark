#!/usr/bin/env python3
"""Which 512-determinant *structure* gives the lowest variational energy?

Compares
  (a) a scattered CIPSI/energy-pruned determinant set,
  (b) a product space  A x B  built from the most important alpha/beta strings
      (the structure SQD / PySCF's selected_ci use), for several aspect ratios,
  (c) a greedily grown product space,
  (d) product space + scattered top-ups.
"""
import os, sys, time, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import Engine, load_ints, key_of, split_key, occ_to_bits
from opt512 import Budget, removal_costs, grow_space

T0 = time.time()


def log(*a):
    print("[%7.1fs]" % (time.time() - T0), *a, flush=True)


def product_keys(astrs, bstrs):
    return np.sort(key_of(np.repeat(astrs, len(bstrs)),
                          np.tile(bstrs, len(astrs))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fcidump")
    ap.add_argument("--grow", type=int, default=80000)
    ap.add_argument("--budget", type=int, default=512)
    ap.add_argument("--outdir", default="out")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    norb, nelec, ms2, h1, g, ecore = load_ints(a.fcidump)
    nocc = (nelec + ms2) // 2
    eng = Engine(norb, nocc, h1, g, ecore)
    B = Budget(eng, a.budget)
    ref = occ_to_bits(range(nocc))
    log("reference determinant E = %.12f" % float(eng.diag([ref], [ref])[0]))

    keys, E, c = grow_space(eng, a.grow, log=log)
    log("CIPSI space |S|=%d  E=%.10f" % (len(keys), E))
    np.savez(os.path.join(a.outdir, "cipsi.npz"), keys=keys, coeff=c, energy=E)

    ka, kb = split_key(keys)
    w = c * c

    # ---------------- (a) scattered ------------------------------------
    o = np.argsort(-np.abs(c))[:6000]
    scat = B.prune(np.sort(keys[o]))
    Es, _ = B.energy(scat)
    log("(a) scattered pruned  %d dets : E = %.12f" % (len(scat), Es))

    # ---------------- (b) product spaces -------------------------------
    ua, inva = np.unique(ka, return_inverse=True)
    ub, invb = np.unique(kb, return_inverse=True)
    wa = np.bincount(inva, weights=w, minlength=len(ua))
    wb = np.bincount(invb, weights=w, minlength=len(ub))
    oa = ua[np.argsort(-wa)]
    ob = ub[np.argsort(-wb)]
    log("distinct alpha strings in CIPSI space: %d" % len(ua))

    best = (Es, scat, "scattered")
    for n in (1, 2, 4, 8, 11, 16, 22, 26, 32, 45, 64, 128, 256, 512):
        m = a.budget // n
        if m < 1 or m > len(ob) or n > len(oa):
            continue
        pk = product_keys(oa[:n], ob[:m])
        Ep, _ = B.energy(pk)
        log("(b) product %3d x %-3d = %4d dets : E = %.12f" % (n, m, len(pk), Ep))
        if Ep < best[0]:
            best = (Ep, pk, "product %dx%d" % (n, m))

    # ---------------- (c) greedy product growth ------------------------
    # start from the best square-ish product and greedily swap strings
    n0 = int(np.sqrt(a.budget))
    A = list(oa[:n0]); Bs = list(ob[:n0])
    pk = product_keys(np.array(A), np.array(Bs))
    Ecur, _ = B.energy(pk)
    log("(c) greedy start %dx%d E=%.12f" % (n0, n0, Ecur))
    pool_a = oa[:400]
    for sweep in range(6):
        improved = False
        for slot in range(len(A)):
            bestg = (Ecur, None)
            for cand in pool_a[:120]:
                if cand in A:
                    continue
                trial = list(A); trial[slot] = cand
                tk = product_keys(np.array(trial), np.array(Bs))
                if len(tk) > a.budget:
                    continue
                Et, _ = B.energy(tk)
                if Et < bestg[0] - 1e-12:
                    bestg = (Et, cand)
            if bestg[1] is not None:
                A[slot] = bestg[1]
                Bs = list(A)          # keep alpha/beta symmetric
                Ecur = bestg[0]
                improved = True
        log("(c) sweep %d: E = %.12f" % (sweep, Ecur))
        if not improved:
            break
    pk = product_keys(np.array(A), np.array(Bs))
    if Ecur < best[0]:
        best = (Ecur, pk, "greedy product")

    # ---------------- (d) product + scattered top-up -------------------
    room = a.budget - len(pk)
    if room > 0:
        Ep, cp = B.energy(pk)
        u, amp, hd, gain = B.candidates(pk, cp, Ep)
        idx = np.argpartition(-gain, room)[:room]
        mix = np.sort(np.concatenate([pk, u[idx]]))
        Em, _ = B.energy(mix)
        log("(d) product+%d scattered = %d dets : E = %.12f" % (room, len(mix), Em))
        if Em < best[0]:
            best = (Em, mix, "product+topup")

    log("BEST structure: %s   E = %.12f" % (best[2], best[0]))
    np.savez(os.path.join(a.outdir, "best_structure.npz"),
             keys=best[1], energy=best[0], label=best[2])


if __name__ == "__main__":
    main()
