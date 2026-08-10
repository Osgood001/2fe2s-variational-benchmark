#!/usr/bin/env python3
"""
Batched Slater-Condon engine on bitmask determinants (alpha_bits, beta_bits).

Determinant phase convention (standard for alpha/beta string CI, = PySCF):
    |D> = (prod_{p in alpha, ascending} a+_{p,alpha})
          (prod_{q in beta , ascending} a+_{q,beta }) |vac>

Everything is vectorised over batches of determinants so that generating the
~7900 single+double connections of tens of thousands of determinants costs
seconds rather than minutes.
"""
import itertools
import numpy as np

NBIT = 20


# ------------------------------------------------------------------ bit utils
_PC16 = np.array([bin(i).count("1") for i in range(1 << 16)], dtype=np.int64)


def popcount(x):
    x = np.asarray(x, dtype=np.int64)
    return _PC16[x & 0xFFFF] + _PC16[(x >> 16) & 0xFFFF]


def occ_to_bits(occ):
    m = 0
    for p in occ:
        m |= 1 << int(p)
    return m


def bitstring(m, norb=NBIT):
    return "".join("1" if (m >> p) & 1 else "0" for p in range(norb))


def bitstring_to_mask(s):
    m = 0
    for i, ch in enumerate(s):
        if ch == "1":
            m |= 1 << i
    return m


def key_of(a, b):
    return (np.asarray(a, dtype=np.int64) << NBIT) | np.asarray(b, dtype=np.int64)


def split_key(k):
    k = np.asarray(k, dtype=np.int64)
    return (k >> NBIT), (k & ((1 << NBIT) - 1))


def load_ints(path):
    from pyscf.tools import fcidump
    from pyscf import ao2mo
    r = fcidump.read(path, verbose=False)
    norb = int(r["NORB"])
    nelec = int(r["NELEC"])
    ms2 = int(r["MS2"])
    h1 = np.asarray(r["H1"], dtype=np.float64)
    g = np.ascontiguousarray(ao2mo.restore(1, np.asarray(r["H2"], np.float64), norb))
    return norb, nelec, ms2, h1, g, float(r["ECORE"])


# ------------------------------------------------------------------ engine
class Engine:
    def __init__(self, norb, nocc, h1, g, ecore):
        self.norb = norb
        self.nocc = nocc
        self.nvir = norb - nocc
        self.h1 = h1
        self.g = g
        self.ecore = ecore
        self.hd = np.diag(h1).copy()
        self.Jd = np.einsum("ppqq->pq", g).copy()
        self.Kd = np.einsum("pqqp->pq", g).copy()
        self.JmK = self.Jd - self.Kd
        self.Jm = np.ascontiguousarray(np.einsum("iapp->iap", g))
        self.Km = np.ascontiguousarray(np.einsum("ippa->iap", g))
        self.JmKm = self.Jm - self.Km

        # string tables ------------------------------------------------
        combs = list(itertools.combinations(range(norb), nocc))
        self.strs = np.array([occ_to_bits(c) for c in combs], dtype=np.int64)
        self.nstr = len(combs)
        self.occ_tab = np.array(combs, dtype=np.int64)
        self.vir_tab = np.array([[p for p in range(norb) if p not in c]
                                 for c in combs], dtype=np.int64)
        self.addr = np.full(1 << norb, -1, dtype=np.int32)
        self.addr[self.strs] = np.arange(self.nstr, dtype=np.int32)
        self.nocc_vec = ((self.strs[:, None] >> np.arange(norb)[None, :]) & 1
                         ).astype(np.float64)
        # same-spin pair index templates
        self.oi, self.oj = np.triu_indices(nocc, 1)
        self.vi, self.vj = np.triu_indices(self.nvir, 1)

    # -------------------------------------------------- diagonal energies
    def diag(self, A, B):
        A = np.atleast_1d(np.asarray(A, dtype=np.int64))
        B = np.atleast_1d(np.asarray(B, dtype=np.int64))
        na = self.nocc_vec[self.addr[A]]
        nb = self.nocc_vec[self.addr[B]]
        e = na @ self.hd + nb @ self.hd
        e += 0.5 * np.einsum("np,pq,nq->n", na, self.JmK, na, optimize=True)
        e += 0.5 * np.einsum("np,pq,nq->n", nb, self.JmK, nb, optimize=True)
        e += np.einsum("np,pq,nq->n", na, self.Jd, nb, optimize=True)
        return e + self.ecore

    def diag_keys(self, keys):
        a, b = split_key(keys)
        return self.diag(a, b)

    # -------------------------------------------------- phase helper
    @staticmethod
    def _sign(mask, i, a):
        lo = np.minimum(i, a)
        hi = np.maximum(i, a)
        one = np.int64(1)
        between = ((one << hi) - 1) & ~((one << (lo + 1)) - 1)
        return 1 - 2 * (popcount(mask & between) & 1)

    # -------------------------------------------------- batched connections
    def connected_batch(self, A, B):
        """For each determinant n in the batch return every determinant d
        reachable by one single or double excitation together with <d|H|n>.

        Returns (keys, vals, src) flat arrays; src is the batch index."""
        A = np.asarray(A, dtype=np.int64)
        B = np.asarray(B, dtype=np.int64)
        n = len(A)
        ia = self.addr[A]
        ib = self.addr[B]
        oa = self.occ_tab[ia]; va = self.vir_tab[ia]      # (n,nocc) (n,nvir)
        ob = self.occ_tab[ib]; vb = self.vir_tab[ib]
        na = self.nocc_vec[ia]; nbv = self.nocc_vec[ib]
        one = np.int64(1)
        no, nv = self.nocc, self.nvir

        keys, vals, srcs = [], [], []
        rng = np.arange(n)

        # ---------- generalised Fock-like intermediates (n,norb,norb)
        Fa = (self.h1[None] + np.einsum("iap,np->nia", self.JmKm, na, optimize=True)
              + np.einsum("iap,np->nia", self.Jm, nbv, optimize=True))
        Fb = (self.h1[None] + np.einsum("iap,np->nia", self.JmKm, nbv, optimize=True)
              + np.einsum("iap,np->nia", self.Jm, na, optimize=True))

        # ---------- alpha singles (n,no,nv)
        I = oa[:, :, None]; Av = va[:, None, :]
        v = Fa[rng[:, None, None], I, Av] * self._sign(A[:, None, None], I, Av)
        nA = A[:, None, None] ^ (one << I) ^ (one << Av)
        keys.append(key_of(nA.ravel(), np.repeat(B, no * nv)))
        vals.append(v.ravel())
        srcs.append(np.repeat(rng, no * nv))

        # ---------- beta singles
        I = ob[:, :, None]; Av = vb[:, None, :]
        v = Fb[rng[:, None, None], I, Av] * self._sign(B[:, None, None], I, Av)
        nB = B[:, None, None] ^ (one << I) ^ (one << Av)
        keys.append(key_of(np.repeat(A, no * nv), nB.ravel()))
        vals.append(v.ravel())
        srcs.append(np.repeat(rng, no * nv))

        # ---------- same-spin doubles
        npo, npv = len(self.oi), len(self.vi)
        for occ, vir, mask, other, alpha in ((oa, va, A, B, True),
                                             (ob, vb, B, A, False)):
            I = occ[:, self.oi][:, :, None]
            J = occ[:, self.oj][:, :, None]
            Av = vir[:, self.vi][:, None, :]
            Bv = vir[:, self.vj][:, None, :]
            v = self.g[I, Av, J, Bv] - self.g[I, Bv, J, Av]
            m0 = mask[:, None, None]
            s1 = self._sign(m0, I, Av)
            m1 = m0 ^ (one << I) ^ (one << Av)
            s2 = self._sign(m1, J, Bv)
            v = v * s1 * s2
            nm = (m1 ^ (one << J) ^ (one << Bv))
            nm = np.broadcast_to(nm, v.shape).ravel()
            if alpha:
                keys.append(key_of(nm, np.repeat(other, npo * npv)))
            else:
                keys.append(key_of(np.repeat(other, npo * npv), nm))
            vals.append(v.ravel())
            srcs.append(np.repeat(rng, npo * npv))

        # ---------- opposite-spin doubles (n,no,nv,no,nv)
        I = oa[:, :, None, None, None]
        Av = va[:, None, :, None, None]
        J = ob[:, None, None, :, None]
        Bv = vb[:, None, None, None, :]
        v = self.g[I, Av, J, Bv]
        sa = self._sign(A[:, None, None, None, None], I, Av)
        sb = self._sign(B[:, None, None, None, None], J, Bv)
        v = v * sa * sb
        nA = np.broadcast_to(A[:, None, None, None, None] ^ (one << I) ^ (one << Av),
                             v.shape).ravel()
        nB = np.broadcast_to(B[:, None, None, None, None] ^ (one << J) ^ (one << Bv),
                             v.shape).ravel()
        keys.append(key_of(nA, nB))
        vals.append(v.ravel())
        srcs.append(np.repeat(rng, (no * nv) ** 2))

        return (np.concatenate(keys), np.concatenate(vals), np.concatenate(srcs))

    # -------------------------------------------------- H|psi> on the outside
    def sigma_external(self, keys, coeff, block=400):
        """(unique_keys, amplitudes) with amplitude_d = sum_I c_I <d|H|I>,
        over every d connected to the space (space members included)."""
        A, B = split_key(np.asarray(keys))
        ks, vs = [], []
        for s in range(0, len(A), block):
            e = min(s + block, len(A))
            k, v, src = self.connected_batch(A[s:e], B[s:e])
            v = v * coeff[s:e][src]
            u, inv = np.unique(k, return_inverse=True)
            ks.append(u)
            vs.append(np.bincount(inv, weights=v, minlength=len(u)))
        k = np.concatenate(ks); v = np.concatenate(vs)
        u, inv = np.unique(k, return_inverse=True)
        return u, np.bincount(inv, weights=v, minlength=len(u))

    # -------------------------------------------------- H inside a space
    def build_H(self, keys, block=400, dense=None):
        import scipy.sparse as sp
        keys = np.asarray(keys, dtype=np.int64)
        n = len(keys)
        order = np.argsort(keys)
        sk = keys[order]
        A, B = split_key(keys)
        rows, cols, data = [], [], []
        for s in range(0, n, block):
            e = min(s + block, n)
            k, v, src = self.connected_batch(A[s:e], B[s:e])
            pos = np.searchsorted(sk, k)
            np.clip(pos, 0, n - 1, out=pos)
            hit = sk[pos] == k
            if hit.any():
                rows.append(src[hit] + s)
                cols.append(order[pos[hit]])
                data.append(v[hit])
        d = self.diag(A, B)
        rows.append(np.arange(n)); cols.append(np.arange(n)); data.append(d)
        rows = np.concatenate(rows); cols = np.concatenate(cols)
        data = np.concatenate(data)
        if dense is None:
            dense = n <= 4000
        if dense:
            H = np.zeros((n, n))
            np.add.at(H, (rows, cols), data)
            return H
        return sp.csr_matrix((data, (rows, cols)), shape=(n, n))

    def solve(self, keys, v0=None, tol=1e-12):
        import scipy.sparse.linalg as spla
        n = len(keys)
        H = self.build_H(keys)
        if isinstance(H, np.ndarray):
            w, V = np.linalg.eigh(H)
            return float(w[0]), V[:, 0]
        kw = dict(k=1, which="SA", tol=tol, maxiter=50000)
        if v0 is not None:
            kw["v0"] = v0
        w, V = spla.eigsh(H, **kw)
        return float(w[0]), V[:, 0]


# ------------------------------------------------------------------ io
def write_csv(path, keys, coeff, norb=NBIT):
    a, b = split_key(np.asarray(keys))
    o = np.argsort(-np.abs(coeff))
    with open(path, "w") as f:
        f.write("alpha,beta,coefficient\n")
        for n in o:
            f.write("%s,%s,%.17g\n" % (bitstring(int(a[n]), norb),
                                       bitstring(int(b[n]), norb),
                                       float(coeff[n])))
