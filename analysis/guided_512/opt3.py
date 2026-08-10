#!/usr/bin/env python3
"""
Strong local search for the lowest-energy B-determinant subspace.

Exact move evaluation, no perturbative approximations in the accept test:

  add   : the energy of S u {d} is the lowest root of the arrow-matrix secular
          equation   H_dd - mu - sum_k b_k^2/(lambda_k - mu) = 0
          with b = V^T <d|H|S>.  Evaluated for *all* candidates at once.
  drop  : the energy of S \ {e} is the lowest root of the Cauchy secular
          equation   sum_k V[e,k]^2/(lambda_k - mu) = 0.
  swap  : add the best d exactly, then drop the exactly cheapest e from the
          resulting B+1 space.  Best-improvement over many candidate d.

Plus k-swaps for coarse progress and large-neighbourhood restarts to escape
local optima.
"""
import os, sys, time, argparse, hashlib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import Engine, load_ints, key_of, split_key, occ_to_bits, write_csv
from opt512 import removal_costs, add_gains, grow_space

T0 = time.time()


def log(*a):
    print("[%8.1fs]" % (time.time() - T0), *a, flush=True)


# ---------------------------------------------------------------- exact add
def add_energies(lam, B, alpha, nbis=90):
    """Lowest eigenvalue of [[diag(lam), b],[b^T, alpha]] for many b rows.

    B is (ncand, n) already rotated into the eigenbasis, alpha is (ncand,).
    The root lies below lam[0]."""
    W = B * B
    n = len(lam)
    lo = np.full(len(alpha), lam[0] - 1.0)
    hi = np.full(len(alpha), lam[0] - 1e-14)

    def f(mu):
        return alpha - mu - (W / (lam[None, :] - mu[:, None])).sum(1)

    for _ in range(60):                      # widen until f(lo) > 0
        bad = f(lo) < 0
        if not bad.any():
            break
        lo = np.where(bad, lam[0] - 2.0 * (lam[0] - lo), lo)
    for _ in range(nbis):
        mid = 0.5 * (lo + hi)
        pos = f(mid) > 0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    return 0.5 * (lo + hi)


class Opt:
    def __init__(self, eng, budget, ntry=32):
        self.eng = eng
        self.B = budget
        self.ntry = ntry
        self.neval = 0

    def solve(self, S):
        H = self.eng.build_H(S, dense=True)
        lam, V = np.linalg.eigh(H)
        self.neval += 1
        return H, lam, V

    def energy(self, S):
        _, lam, V = self.solve(S)
        return float(lam[0]), V[:, 0]

    # ---- couplings <d|H|I> for all candidates d and all I in S
    def couplings(self, S, cand):
        A, Bt = split_key(S)
        order = np.argsort(cand)
        sc = cand[order]
        C = np.zeros((len(cand), len(S)))
        for s in range(0, len(S), 400):
            e = min(s + 400, len(S))
            k, v, src = self.eng.connected_batch(A[s:e], Bt[s:e])
            pos = np.searchsorted(sc, k)
            np.clip(pos, 0, len(sc) - 1, out=pos)
            hit = sc[pos] == k
            np.add.at(C, (order[pos[hit]], src[hit] + s), v[hit])
        return C

    def gen_candidates(self, S, c, E, cap=6000):
        u, amp = self.eng.sigma_external(S, c)
        m = ~np.isin(u, S)
        u, amp = u[m], amp[m]
        hd = self.eng.diag_keys(u)
        gain = add_gains(E, hd, amp)
        if len(u) > cap:
            i = np.argpartition(-gain, cap)[:cap]
            u, hd, gain = u[i], hd[i], gain[i]
        o = np.argsort(-gain)
        return u[o], hd[o]

    # ---- one best-improvement exact swap sweep -------------------------
    def swap_sweep(self, S, ntry=None, cap=6000, deadline=None):
        ntry = ntry or self.ntry
        H, lam, V = self.solve(S)
        E = float(lam[0])
        cand, hdd = self.gen_candidates(S, V[:, 0], E, cap=cap)
        if len(cand) == 0:
            return S, E, False
        C = self.couplings(S, cand)
        Brot = C @ V                       # (ncand, B) in eigenbasis
        eadd = add_energies(lam, Brot, hdd)
        order = np.argsort(eadd)
        best = (E, None, None)
        for j in order[:ntry]:
            if deadline and time.time() > deadline:
                break
            n = self.B
            Hp = np.empty((n + 1, n + 1))
            Hp[:n, :n] = H
            Hp[:n, n] = C[j]
            Hp[n, :n] = C[j]
            Hp[n, n] = hdd[j]
            lam2, V2 = np.linalg.eigh(Hp)
            self.neval += 1
            newE = removal_costs(lam2, V2)
            newE[n] = np.inf                      # do not remove what we added
            e = int(np.argmin(newE))
            if newE[e] < best[0] - 1e-13:
                best = (float(newE[e]), int(j), e)
        if best[1] is None:
            return S, E, False
        j, e = best[1], best[2]
        S2 = np.sort(np.concatenate([np.delete(S, e), [cand[j]]]))
        return S2, best[0], True

    # ---- coarse k-swap -------------------------------------------------
    def kswap(self, S, deadline=None):
        H, lam, V = self.solve(S)
        E = float(lam[0])
        cost = removal_costs(lam, V) - E
        cand, hdd = self.gen_candidates(S, V[:, 0], E)
        oc = np.argsort(cost)
        for k in (128, 64, 32, 16, 8, 4, 2):
            if k > len(cand):
                continue
            trial = np.sort(np.concatenate([np.delete(S, oc[:k]), cand[:k]]))
            if len(np.unique(trial)) != self.B:
                continue
            E2, _ = self.energy(trial)
            if E2 < E - 1e-13:
                return trial, E2, True
        return S, E, False

    def descend(self, S, deadline, tag="", report=25):
        E, _ = self.energy(S)
        it = 0
        while True:
            if deadline and time.time() > deadline:
                break
            S2, E2, ok = self.kswap(S, deadline)
            if ok and E2 < E - 1e-13:
                S, E = S2, E2
            else:
                S2, E2, ok = self.swap_sweep(S, deadline=deadline)
                if not ok or E2 >= E - 1e-13:
                    break
                S, E = S2, E2
            it += 1
            if it % report == 0:
                log("   %s it=%4d E=%.12f" % (tag, it, E))
        return S, E


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fcidump")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--budget", type=int, default=512)
    ap.add_argument("--pool", default=None)
    ap.add_argument("--start", default=None)
    ap.add_argument("--grow", type=int, default=120000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--hours", type=float, default=2.0)
    ap.add_argument("--restarts", type=int, default=1000)
    ap.add_argument("--ntry", type=int, default=32)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    deadline = time.time() + a.hours * 3600
    rng = np.random.default_rng(a.seed)

    log("fcidump sha256", hashlib.sha256(open(a.fcidump, "rb").read()).hexdigest())
    norb, nelec, ms2, h1, g, ecore = load_ints(a.fcidump)
    nocc = (nelec + ms2) // 2
    eng = Engine(norb, nocc, h1, g, ecore)
    ref = occ_to_bits(range(nocc))
    e_ref = float(eng.diag([ref], [ref])[0])
    log("NORB=%d NELEC=%d MS2=%d  closed-shell ref E=%.12f" % (norb, nelec, ms2, e_ref))
    O = Opt(eng, a.budget, ntry=a.ntry)

    if a.pool and os.path.exists(a.pool):
        d = np.load(a.pool, allow_pickle=True)
        if "strings" in d:
            strs = d["strings"].astype(np.int64)
            pk = key_of(strs[d["alpha_addr"]], strs[d["beta_addr"]])
            pc = np.asarray(d["coeff"])
            log("FCI pool %d dets  E_FCI=%.12f" % (len(pk), float(d["e_fci"])))
        else:
            pk, pc = np.asarray(d["keys"]), np.asarray(d["coeff"])
            log("pool %d dets" % len(pk))
    else:
        pk, Eg, pc = grow_space(eng, a.grow, log=log)
        np.savez(os.path.join(a.outdir, "cipsi.npz"), keys=pk, coeff=pc, energy=Eg)
    o = np.argsort(-np.abs(pc))
    pk, pc = pk[o], pc[o]

    if a.start and os.path.exists(a.start):
        S = np.sort(np.load(a.start)["keys"])
        log("warm start %s" % a.start)
    else:
        S = np.sort(pk[:a.budget])
    E, _ = O.energy(S)
    log("start E = %.12f" % E)

    S, E = O.descend(S, deadline, tag="D0")
    log("after descent: E = %.12f" % E)
    best_S, best_E = S.copy(), E
    np.savez(os.path.join(a.outdir, "best.npz"), keys=best_S, energy=best_E)

    # ---------------- large-neighbourhood restarts --------------------
    for r in range(a.restarts):
        if time.time() > deadline:
            break
        _, lam, V = O.solve(best_S)
        c = V[:, 0]
        k = int(rng.integers(3, 48))
        p = 1.0 / (np.abs(c) + 1e-8)
        p /= p.sum()
        drop = rng.choice(len(best_S), size=k, replace=False, p=p)
        partial = np.delete(best_S, drop)
        Ep, cp = O.energy(partial)
        cand, hdd = O.gen_candidates(partial, cp, Ep)
        pool_extra = pk[:400000]
        pool_extra = pool_extra[~np.isin(pool_extra, partial)]
        if len(pool_extra) and rng.random() < 0.5:
            take_pool = min(len(pool_extra), 4000)
            sel = rng.choice(len(pool_extra), size=take_pool, replace=False)
            cand = np.concatenate([cand, pool_extra[sel]])
        cand = cand[~np.isin(cand, partial)]
        need = a.budget - len(partial)
        pick = rng.choice(min(len(cand), max(6 * need, 120)), size=need, replace=False)
        S = np.sort(np.concatenate([partial, cand[pick]]))
        if len(np.unique(S)) != a.budget:
            continue
        S, E = O.descend(S, deadline, tag="R%d" % r, report=10**9)
        if E < best_E - 1e-13:
            best_S, best_E = S.copy(), E
            log("restart %4d k=%-3d ACCEPT  E=%.12f" % (r, k, best_E))
            np.savez(os.path.join(a.outdir, "best.npz"), keys=best_S, energy=best_E)
        elif r % 20 == 0:
            log("restart %4d k=%-3d  (%.3e above best)" % (r, k, E - best_E))

    Ef, cf = O.energy(best_S)
    log("FINAL %d-determinant energy: %.12f Eh  (%.6f below closed shell)  evals=%d"
        % (a.budget, Ef, Ef - e_ref, O.neval))
    np.savez(os.path.join(a.outdir, "best.npz"), keys=best_S, coeff=cf, energy=Ef)
    write_csv(os.path.join(a.outdir, "wavefunction.csv"), best_S, cf, norb)
    log("wrote", os.path.join(a.outdir, "wavefunction.csv"))


if __name__ == "__main__":
    main()
