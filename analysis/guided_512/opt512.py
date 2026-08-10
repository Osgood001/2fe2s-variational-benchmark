#!/usr/bin/env python3
"""
Fixed-budget variational optimiser: find the 512-determinant subspace whose
lowest eigenvalue of H is as low as possible, and report its eigenvector.

Ingredients
-----------
* CIPSI/heat-bath growth to build a good starting space and candidate pool.
* Exact cost of deleting one determinant from a space, via the Cauchy secular
  equation  sum_k V[e,k]^2 / (lambda_k - mu) = 0  (no extra diagonalisation).
* Exact gain of appending one determinant d to the current psi (2x2 problem,
  |d> is orthogonal to psi because d is outside the space).
* grow -> prune cycles plus large-neighbourhood-search swaps.
"""
import os, sys, time, argparse, hashlib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from engine import (Engine, load_ints, key_of, split_key, occ_to_bits,
                    bitstring, write_csv)

T0 = time.time()


def log(*a):
    print("[%8.1fs]" % (time.time() - T0), *a, flush=True)


# ------------------------------------------------------------------ exact tools
def removal_costs(lam, V, eps=1e-13):
    """Exact  min_{x_e=0} <x|H|x>  for every e, given H = V diag(lam) V^T.

    Returns the array of new ground-state energies (one per deleted index)."""
    n = len(lam)
    W = V ** 2                     # W[e,k] = V[e,k]^2
    lo = np.full(n, lam[0])
    hi = np.full(n, lam[1] if n > 1 else lam[0] + 1.0)
    # guard against exact degeneracy of the two lowest roots
    if n > 1 and lam[1] - lam[0] < 1e-12:
        j = 1
        while j < n and lam[j] - lam[0] < 1e-12:
            j += 1
        hi = np.full(n, lam[min(j, n - 1)])
    lo = lo + eps
    hi = hi - eps
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = (W / (lam[None, :] - mid[:, None])).sum(1)
        neg = f < 0.0
        lo = np.where(neg, mid, lo)
        hi = np.where(neg, hi, mid)
    return 0.5 * (lo + hi)


def add_gains(E, hd, amp):
    """Exact lowering from appending one external determinant to fixed psi."""
    d = hd - E
    return 0.5 * (np.sqrt(d * d + 4.0 * amp * amp) - d)


class Budget:
    def __init__(self, eng, budget):
        self.eng = eng
        self.B = budget
        self.n_solve = 0

    def dense_solve(self, keys):
        H = self.eng.build_H(keys, dense=True)
        lam, V = np.linalg.eigh(H)
        self.n_solve += 1
        return lam, V

    def energy(self, keys):
        lam, V = self.dense_solve(keys)
        return float(lam[0]), V[:, 0]

    # ---- candidate generation --------------------------------------
    def candidates(self, keys, c, E, nmax=None):
        u, amp = self.eng.sigma_external(keys, c)
        inside = np.isin(u, keys)
        u = u[~inside]; amp = amp[~inside]
        hd = self.eng.diag_keys(u)
        gain = add_gains(E, hd, amp)
        if nmax is not None and len(u) > nmax:
            idx = np.argpartition(-gain, nmax)[:nmax]
            u, amp, hd, gain = u[idx], amp[idx], hd[idx], gain[idx]
        return u, amp, hd, gain

    # ---- prune a space down to the budget --------------------------
    def prune(self, keys):
        keys = np.asarray(keys)
        while len(keys) > self.B:
            lam, V = self.dense_solve(keys)
            newE = removal_costs(lam, V)
            excess = len(keys) - self.B
            drop = max(1, min(excess, int(np.ceil(0.30 * len(keys)))))
            drop = min(drop, excess) if excess < 0.30 * len(keys) else drop
            order = np.argsort(newE)             # cheapest to delete first
            keep = np.ones(len(keys), dtype=bool)
            keep[order[:drop]] = False
            keys = keys[keep]
        return keys


def grow_space(eng, target, log=log, start=None):
    """CIPSI growth from the closed-shell determinant (or a given start)."""
    if start is None:
        ref = occ_to_bits(range(eng.nocc))
        keys = key_of([ref], [ref])
        c = np.array([1.0])
        E = float(eng.diag([ref], [ref])[0])
    else:
        keys = start
        E, c = eng.solve(keys)
    while len(keys) < target:
        u, amp = eng.sigma_external(keys, c)
        mask = ~np.isin(u, keys)
        u, amp = u[mask], amp[mask]
        if len(u) == 0:
            break
        hd = eng.diag_keys(u)
        den = E - hd
        den[np.abs(den) < 1e-9] = -1e-9
        e2 = amp * amp / den
        ntake = min(len(u), max(4000, int(1.6 * len(keys))), target - len(keys) + 1)
        idx = np.argpartition(e2, ntake - 1)[:ntake]
        keys = np.sort(np.concatenate([keys, u[idx]]))
        E, c = eng.solve(keys)
        log("   grow |S|=%-8d E=%.10f   PT2(discarded)=%.6f" %
            (len(keys), E, float(e2.sum() - e2[idx].sum())))
    return keys, E, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fcidump")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--budget", type=int, default=512)
    ap.add_argument("--grow", type=int, default=150000)
    ap.add_argument("--pool", default=None)
    ap.add_argument("--start", default=None, help="npz with keys= for a warm start")
    ap.add_argument("--cycles", type=int, default=40)
    ap.add_argument("--add", type=int, default=2048)
    ap.add_argument("--lns", type=int, default=0, help="large-neighbourhood rounds")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    log("fcidump sha256",
        hashlib.sha256(open(args.fcidump, "rb").read()).hexdigest())
    norb, nelec, ms2, h1, g, ecore = load_ints(args.fcidump)
    nocc = (nelec + ms2) // 2
    log("NORB=%d NELEC=%d MS2=%d -> (%d,%d)  ECORE=%.12f" %
        (norb, nelec, ms2, nocc, (nelec - ms2) // 2, ecore))
    eng = Engine(norb, nocc, h1, g, ecore)
    ref = occ_to_bits(range(nocc))
    e_ref = float(eng.diag([ref], [ref])[0])
    log("closed-shell reference determinant E = %.12f" % e_ref)

    B = Budget(eng, args.budget)

    # -------------------------------------------------- starting space
    if args.start and os.path.exists(args.start):
        d = np.load(args.start)
        cur = np.sort(d["keys"])
        log("warm start from %s: %d determinants" % (args.start, len(cur)))
    else:
        if args.pool and os.path.exists(args.pool):
            d = np.load(args.pool, allow_pickle=True)
            if "strings" in d:
                strs = d["strings"].astype(np.int64)
                pk = key_of(strs[d["alpha_addr"]], strs[d["beta_addr"]])
                pc = d["coeff"]
            else:
                pk, pc = d["keys"], d["coeff"]
            o = np.argsort(-np.abs(pc))
            big = np.sort(pk[o[:max(args.budget * 40, 20000)]])
            log("pool %d determinants -> working space %d" % (len(pk), len(big)))
            cur = B.prune(big)
        else:
            keys, E, c = grow_space(eng, args.grow)
            log("grown space |S|=%d  E=%.12f" % (len(keys), E))
            np.savez(os.path.join(args.outdir, "grown.npz"), keys=keys, coeff=c,
                     energy=E)
            o = np.argsort(-np.abs(c))
            cur = B.prune(np.sort(keys[o[:min(len(keys), args.budget * 12)]]))

    cur = np.sort(cur)
    E, c = B.energy(cur)
    log("start: %d determinants, E = %.12f" % (len(cur), E))
    best_E, best_keys, best_c = E, cur.copy(), c.copy()

    # -------------------------------------------------- grow / prune cycles
    stall = 0
    for it in range(args.cycles):
        u, amp, hd, gain = B.candidates(cur, c, E)
        m = min(args.add, len(u))
        idx = np.argpartition(-gain, m - 1)[:m]
        space = np.unique(np.concatenate([cur, u[idx]]))
        Ebig, _ = eng.solve(space)
        new = B.prune(space)
        Enew, cnew = B.energy(new)
        log("cycle %2d |P|=%5d E(P)=%.10f -> |S|=%d E=%.12f  (dE=%+.3e)"
            % (it, len(space), Ebig, len(new), Enew, Enew - E))
        if Enew < best_E - 1e-13:
            best_E, best_keys, best_c = Enew, new.copy(), cnew.copy()
        if Enew > E - 1e-11:
            stall += 1
        else:
            stall = 0
        cur, E, c = new, Enew, cnew
        if stall >= 2:
            log("grow/prune converged")
            break

    # -------------------------------------------------- LNS polish
    cur, E, c = best_keys.copy(), best_E, best_c.copy()
    for it in range(args.lns):
        lam, V = B.dense_solve(cur)
        E = float(lam[0]); c = V[:, 0]
        newE = removal_costs(lam, V)
        u, amp, hd, gain = B.candidates(cur, c, E)
        k = int(rng.integers(8, 96))
        weak = np.argsort(newE)[:k]                     # cheapest to remove
        strong = np.argpartition(-gain, k)[:k]
        trial = np.unique(np.concatenate([np.delete(cur, weak), u[strong]]))
        if len(trial) > args.budget:
            trial = B.prune(trial)
        Et, ct = B.energy(trial)
        if Et < E - 1e-13:
            cur, E, c = trial, Et, ct
            if Et < best_E:
                best_E, best_keys, best_c = Et, trial.copy(), ct.copy()
            log("  lns %3d k=%-3d accept E=%.12f" % (it, k, Et))
        elif it % 10 == 0:
            log("  lns %3d k=%-3d reject (%.3e)" % (it, k, Et - E))

    log("BEST %d-determinant variational energy: %.12f Eh" % (args.budget, best_E))
    log("   vs closed-shell reference: %.6f Eh" % (best_E - e_ref))
    np.savez(os.path.join(args.outdir, "best.npz"),
             keys=best_keys, coeff=best_c, energy=best_E)
    write_csv(os.path.join(args.outdir, "wavefunction.csv"), best_keys, best_c, norb)
    log("wrote", os.path.join(args.outdir, "wavefunction.csv"))


if __name__ == "__main__":
    main()
