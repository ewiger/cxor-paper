"""
Verification harness for the cxor paper (doc/prompt.md spec).

Notation follows the spec:
    x   tape (input),  x: K -> F_2
    s   table (state), s: A -> F_2,  A = F_2^m,  |A| = 2^m
    pi  prefix (initial address), pi in F_2^m
    w   wiring (write rule), w: F_2 x F_2 -> F_2
    r   output
    cxor^w_m(s, pi, x) = r

Wiring index:  W_k,  k = 8*w(0,0) + 4*w(0,1) + 2*w(1,0) + 1*w(1,1)
ANF:           w(s,y) = d0 + d1*s + d2*y + d3*s*y
"""
from itertools import product
import random

# ------------------------------------------------------------------ wirings

def w_tt(k):
    """Truth table [w(0,0), w(0,1), w(1,0), w(1,1)] of W_k."""
    return [(k >> (3 - i)) & 1 for i in range(4)]


def w_of(k):
    t = w_tt(k)
    return lambda s, y: t[2 * s + y]


def anf(k):
    """(d0,d1,d2,d3) with w = d0 + d1 s + d2 y + d3 s y  (Mobius inversion)."""
    w = w_of(k)
    d0 = w(0, 0)
    d1 = w(1, 0) ^ d0
    d2 = w(0, 1) ^ d0
    d3 = w(1, 1) ^ d0 ^ d1 ^ d2
    return d0, d1, d2, d3


def anf_str(k):
    t = [n for n, on in zip(["1", "s", "y", "sy"], anf(k)) if on]
    return " + ".join(t) if t else "0"


def index_of(w):
    return 8 * w(0, 0) + 4 * w(0, 1) + 2 * w(1, 0) + 1 * w(1, 1)


# ------------------------------------------------------------------ the machine

def cxor(s, pi, x, k=5, m=None):
    """cxor^w_m(s, pi, x) -> (r, addresses).  Window EXCLUDES the current symbol."""
    if m is None:
        m = len(pi)
    w = w_of(k)
    s = list(s)
    cxt = list(pi)
    r, addr = [], []
    for sym in x:
        a = int("".join(map(str, cxt)), 2) if m else 0
        addr.append(a)
        cell = s[a]
        r.append(sym ^ cell)
        s[a] = w(cell, sym)
        if m:
            cxt.pop(0); cxt.append(sym)
    return r, addr


def decode(s, pi, r, k=5, m=None):
    """Directional inverse: regenerates addressing from recovered symbols."""
    if m is None:
        m = len(pi)
    w = w_of(k)
    s = list(s)
    cxt = list(pi)
    out = []
    for rk in r:
        a = int("".join(map(str, cxt)), 2) if m else 0
        cell = s[a]
        sym = rk ^ cell
        out.append(sym)
        s[a] = w(cell, sym)
        if m:
            cxt.pop(0); cxt.append(sym)
    return out


def frozen(s, addr, N, x, k=5):
    """Run with the address sequence supplied as external data."""
    w = w_of(k)
    s = list(s)
    r = []
    for sym, a in zip(x, addr):
        cell = s[a]
        r.append(sym ^ cell)
        s[a] = w(cell, sym)
    return r


# ------------------------------------------------------------------ Post's clones

def in_T0(k): return w_of(k)(0, 0) == 0
def in_T1(k): return w_of(k)(1, 1) == 1
def in_A(k):  return anf(k)[3] == 0


def in_M(k):
    w = w_of(k)
    pts = list(product([0, 1], repeat=2))
    return all(w(*p) <= w(*q) for p in pts for q in pts
               if p[0] <= q[0] and p[1] <= q[1])


def in_D(k):
    w = w_of(k)
    return all(w(1 ^ a, 1 ^ b) == 1 ^ w(a, b) for a in (0, 1) for b in (0, 1))


CLONES = [("T0", in_T0), ("T1", in_T1), ("M", in_M), ("D", in_D), ("A", in_A)]


def clones_of(k):
    return [nm for nm, f in CLONES if f(k)]


# ------------------------------------------------------------------ F_2[[T]]

def series_apply(num, den, x, prev):
    """Apply the rational operator num(T)/den(T) to x, with T the predecessor shift
    given by prev[k] (index of the previous position in k's class, or None).

    Solves den(T) r = num(T) x coordinatewise; den must have constant term 1."""
    n = len(x)
    r = [0] * n
    for k in range(n):
        acc = 0
        for j, c in enumerate(num):          # num[j] is the coefficient of T^j
            if c:
                idx = k
                for _ in range(j):
                    idx = prev[idx] if idx is not None else None
                    if idx is None:
                        break
                if idx is not None:
                    acc ^= x[idx]
        for j, c in enumerate(den):
            if c and j > 0:
                idx = k
                for _ in range(j):
                    idx = prev[idx] if idx is not None else None
                    if idx is None:
                        break
                if idx is not None:
                    acc ^= r[idx]
        r[k] = acc
    return r


def prev_map(addr):
    """prev[k] = previous position in k's address class, or None."""
    last, prev = {}, []
    for k, a in enumerate(addr):
        prev.append(last.get(a))
        last[a] = k
    return prev


# =============================================================================
# Checks.  Each returns nothing and prints; `python3 cta_spec.py` runs them all.
# Section numbers refer to the paper.
# =============================================================================

def check_index_and_lattice():
    """S4 -- index convention, ANF, Post's clones, functional completeness."""
    print("wiring index W_k, k = 8w(0,0)+4w(0,1)+2w(1,0)+1w(1,1)")
    ok = all(index_of(w_of(k)) == k for k in range(16))
    print(f"  index round-trips for all 16: {ok}")
    print(f"  W5 = {anf_str(5)!r} (projection), clones {clones_of(5)}")
    print(f"  W6 = {anf_str(6)!r} (XOR),        clones {clones_of(6)}")
    A = [k for k in range(16) if in_A(k)]
    print(f"  affine wirings (d3=0): {A}")
    print(f"  W8 is NOR: {[w_of(8)(a,b) for a in (0,1) for b in (0,1)] == [1,0,0,0]}; "
          f"W14 is NAND: {[w_of(14)(a,b) for a in (0,1) for b in (0,1)] == [1,1,1,0]}")
    print(f"  in no maximal clone (functionally complete alone): "
          f"{[k for k in range(16) if not clones_of(k)]}")


def check_causality_and_bijectivity():
    """S5 -- r(k) depends only on x(0..k); every wiring is bijective in the tape."""
    bad = []
    for k in range(16):
        for _ in range(200):
            m = random.choice([0, 1, 2, 3]); N = 1 << m; L = random.randint(2, 40)
            pi = [random.randint(0, 1) for _ in range(m)]
            s = [random.randint(0, 1) for _ in range(N)]
            x = [random.randint(0, 1) for _ in range(L)]
            j = random.randrange(L)
            y = list(x); y[j] ^= 1
            if cxor(s, pi, x, k, m)[0][:j] != cxor(s, pi, y, k, m)[0][:j]:
                bad.append(("causality", k)); break
            if decode(s, pi, cxor(s, pi, x, k, m)[0], k, m) != x:
                bad.append(("decode", k)); break
    print(f"  causality and decode(cxor(x)) = x, all 16 wirings: {not bad}"
          + (f"  failures {bad}" if bad else ""))

    for m, L in ((0, 12), (1, 10), (2, 8)):
        N = 1 << m
        bad = []
        for k in range(16):
            for s in product([0, 1], repeat=N):
                for pi in product([0, 1], repeat=m):
                    p = list(pi)
                    for Lx in range(1, L + 1):
                        for t in product([0, 1], repeat=Lx):
                            x = list(t)
                            if decode(s, p, cxor(s, p, x, k, m)[0], k, m) != x:
                                bad.append((k, s, pi, x)); break
        print(f"    exhaustive at m={m}: all 16 wirings, all {2 ** N} tables, "
              f"all {2 ** m} prefixes, all tapes to length {L}: {not bad}"
              + (f"  failures {bad[:3]}" if bad else ""))
    bad = []
    for k in range(16):
        for m in (1, 2):
            N = 1 << m
            for _ in range(4):
                s = [random.randint(0, 1) for _ in range(N)]
                pi = [random.randint(0, 1) for _ in range(m)]
                for L in range(1, 9):
                    if len({tuple(cxor(s, pi, list(t), k, m)[0])
                            for t in product([0, 1], repeat=L)}) != 2 ** L:
                        bad.append(k)
    print(f"  bijective in the tape, exhaustive, all 16 wirings: {not bad}")

    # The other alignment (window shifted BEFORE the lookup) destroys injectivity.
    def lagging(s, pi, x, k, m):
        w = w_of(k); s = list(s); cxt = list(pi); r = []
        for sym in x:
            if m:
                cxt.pop(0); cxt.append(sym)
            a = int("".join(map(str, cxt)), 2) if m else 0
            cell = s[a]; r.append(sym ^ cell); s[a] = w(cell, sym)
        return r

    coll = [k for k in range(16)
            if lagging([0, 1], [0], [0], k, 1) == lagging([0, 1], [0], [1], k, 1)]
    print(f"  leading window, s=(0,1), m=1: tapes '0' and '1' collide, all 16: "
          f"{len(coll) == 16}")
    surv = []
    for k in range(16):
        seen = {}
        hit = False
        for L in range(1, 7):
            for t in product([0, 1], repeat=L):
                key = tuple(lagging([0, 0], [0], list(t), k, 1))
                if key in seen and seen[key] != t:
                    hit = True; break
                seen[key] = t
            if hit: break
        if not hit:
            surv.append(k)
    print(f"  leading window, zero table: non-injective for {16 - len(surv)}/16; "
          f"survivors {['W%d' % k for k in surv]} (d0 = d2 = 0, the identity map)")


def check_affine_in_table():
    """S6 -- the main theorem, plus confinement, rank and support."""
    bad = []
    for k in range(16):
        for _ in range(300):
            m = random.choice([0, 1, 2, 3]); N = 1 << m; L = random.randint(1, 60)
            pi = [random.randint(0, 1) for _ in range(m)]
            x = [random.randint(0, 1) for _ in range(L)]
            s = [random.randint(0, 1) for _ in range(N)]
            t = [random.randint(0, 1) for _ in range(N)]
            lhs = cxor([p ^ q for p, q in zip(s, t)], pi, x, k, m)[0]
            rhs = [a ^ b ^ c for a, b, c in zip(cxor(s, pi, x, k, m)[0],
                                                cxor(t, pi, x, k, m)[0],
                                                cxor([0]*N, pi, x, k, m)[0])]
            if lhs != rhs:
                bad.append(k); break
    print(f"  cxor(s+s') = cxor(s)+cxor(s')+cxor(0), all 16 wirings: {not bad}")

    bad = tot = 0
    for _ in range(300):
        k = random.randrange(16); m = random.choice([1, 2, 3]); N = 1 << m
        L = random.randint(5, 80)
        pi = [random.randint(0, 1) for _ in range(m)]
        x = [random.randint(0, 1) for _ in range(L)]
        s = [random.randint(0, 1) for _ in range(N)]
        r0, addr = cxor(s, pi, x, k, m)
        for a in range(N):
            s2 = list(s); s2[a] ^= 1
            r1 = cxor(s2, pi, x, k, m)[0]
            tot += 1
            bad += any(r0[j] != r1[j] and addr[j] != a for j in range(L))
    print(f"  confinement: perturbing one cell changes nothing outside its class: "
          f"{bad == 0} ({bad}/{tot})")

    bad = []
    for k in range(16):
        for _ in range(80):
            m = random.choice([1, 2]); N = 1 << m; L = random.randint(1, 30)
            pi = [random.randint(0, 1) for _ in range(m)]
            x = [random.randint(0, 1) for _ in range(L)]
            base, addr = cxor([0]*N, pi, x, k, m)
            rank = 0
            for a in range(N):
                e = [0]*N; e[a] = 1
                if any(p ^ q for p, q in zip(cxor(e, pi, x, k, m)[0], base)):
                    rank += 1
            V = len(set(addr))
            inj = len({tuple(cxor(list(t), pi, x, k, m)[0])
                       for t in product([0, 1], repeat=N)}) == 2 ** N
            if rank != V or inj != (V == N):
                bad.append(k); break
    print(f"  rank P = |V|, and injective in s iff every address is visited: {not bad}")


def check_support():
    """S6 -- supp P is the first-visit set for W5 and all of K for W6."""
    m, L = 4, 400; N = 1 << m
    pi = [0]*m; x = [random.randint(0, 1) for _ in range(L)]
    for k in (5, 6):
        base, addr = cxor([0]*N, pi, x, k, m)
        sup = set()
        for a in range(N):
            e = [0]*N; e[a] = 1
            sup |= {j for j, (p, q) in enumerate(zip(cxor(e, pi, x, k, m)[0], base)) if p != q}
        seen = set(); first = set()
        for j, a in enumerate(addr):
            if a not in seen:
                seen.add(a); first.add(j)
        tag = "= first visits" if sup == first else ("= all of K" if len(sup) == L else "?")
        print(f"  W{k}: |supp P| = {len(sup)}/{L}   {tag}")


def check_not_affine_in_prefix():
    """S6 -- pi changes the addressing, so cxor is not affine in the prefix."""
    for m in (1, 2):
        N = 1 << m
        for L in range(1, 7):
            for x in product([0, 1], repeat=L):
                for p, q in product(product([0, 1], repeat=m), repeat=2):
                    a = cxor([0]*N, list(p), list(x), 5, m)[0]
                    b = cxor([0]*N, list(q), list(x), 5, m)[0]
                    z = cxor([0]*N, [0]*m, list(x), 5, m)[0]
                    d = cxor([0]*N, [u ^ v for u, v in zip(p, q)], list(x), 5, m)[0]
                    if d != [i ^ j ^ l for i, j, l in zip(a, b, z)]:
                        print(f"  counterexample (W5): m={m}, x={''.join(map(str,x))}, "
                              f"pi={p} vs {q}")
                        return


def check_transfer():
    """S7 -- per-class transfer function and the closure theorem."""
    L = 14; PREV = [None] + list(range(L - 1))
    bad = []
    for k in [k for k in range(16) if in_A(k)]:
        d0, d1, d2, _ = anf(k)
        num, den = [1, d1 ^ d2], [1, d1]
        for s0 in (0, 1):
            c = cxor([s0], [], [0]*L, k, 0)[0]
            for _ in range(300):
                x = [random.randint(0, 1) for _ in range(L)]
                want = [a ^ b for a, b in zip(series_apply(num, den, x, PREV), c)]
                if cxor([s0], [], x, k, 0)[0] != want:
                    bad.append(k); break
    print(f"  r = (1+(d1+d2)T)/(1+d1 T) x + c, all 8 affine wirings: {not bad}")
    groups = {}
    for k in [k for k in range(16) if in_A(k)]:
        _, d1, d2, _ = anf(k)
        op = "1" if (d1, d2) in ((0, 0), (1, 0)) else ("1+T" if (d1, d2) == (0, 1)
                                                       else "(1+T)^-1")
        groups.setdefault(op, []).append(k)
    for op in ("1", "1+T", "(1+T)^-1"):
        free = [k for k in groups[op] if anf(k)[0] == 0]
        print(f"    {op:<10} <- W{groups[op]}   constant-free: W{free}")
    print(f"  exactly three operators arise: {len(groups) == 3}")
    nl = []
    for k in [k for k in range(16) if not in_A(k)]:
        lin = True
        for s0 in (0, 1):
            base = cxor([s0], [], [0]*L, k, 0)[0]
            for _ in range(200):
                a = [random.randint(0, 1) for _ in range(L)]
                b = [random.randint(0, 1) for _ in range(L)]
                fa = cxor([s0], [], a, k, 0)[0]; fb = cxor([s0], [], b, k, 0)[0]
                fab = cxor([s0], [], [p ^ q for p, q in zip(a, b)], k, 0)[0]
                if fab != [p ^ q ^ r for p, q, r in zip(fa, fb, base)]:
                    lin = False; break
        if not lin:
            nl.append(k)
    print(f"  the 8 non-affine wirings realise no transfer function (nonlinear in x): "
          f"{nl == [k for k in range(16) if not in_A(k)]}")


def check_gray():
    """S8 -- m=0 is exactly the Gray code; m=1 is stride 2; m>=2 is not linear."""
    L = 12
    b = lambda n: [(n >> (L - 1 - i)) & 1 for i in range(L)]
    v = lambda t: int("".join(map(str, t)), 2)
    ok5 = all(v(cxor([0], [], b(n), 5, 0)[0]) == (n ^ (n >> 1)) for n in range(1 << L))
    ok6 = all(v(cxor([0], [], b(n ^ (n >> 1)), 6, 0)[0]) == n for n in range(1 << L))
    print(f"  m=0: W5 is binary->Gray {ok5}; W6 is Gray->binary {ok6} (all {1<<L} words)")
    ok = True
    for L2 in (1, 2, 3, 8, 13):
        for t in product([0, 1], repeat=L2):
            x = list(t)
            if cxor([0, 0], [0], x, 5, 1)[0] != [x[j] ^ (x[j-2] if j >= 2 else 0)
                                                 for j in range(L2)]:
                ok = False
    print(f"  m=1, s=0, pi=0: W5 collapses to r(k) = x(k) + x(k-2): {ok}")
    for m in (0, 1, 2, 3):
        N = 1 << m; L3 = 12

        def islin(seed):
            for _ in range(200):
                s0 = seed()
                base = cxor(s0, [0]*m, [0]*L3, 5, m)[0]
                a = [random.randint(0, 1) for _ in range(L3)]
                bb = [random.randint(0, 1) for _ in range(L3)]
                fa = cxor(s0, [0]*m, a, 5, m)[0]; fb = cxor(s0, [0]*m, bb, 5, m)[0]
                fab = cxor(s0, [0]*m, [p ^ q for p, q in zip(a, bb)], 5, m)[0]
                if fab != [p ^ q ^ r for p, q, r in zip(fa, fb, base)]:
                    return False
            return True
        print(f"    m={m}: linear in x with s=0 -> {islin(lambda: [0]*N)};   "
              f"with random s -> {islin(lambda: [random.randint(0,1) for _ in range(N)])}")


def check_duality_and_inverses():
    """S9 -- the three notions of inverse, and the visit-index bound."""
    bad = 0
    for _ in range(2000):
        N = random.randint(1, 10); L = random.randint(1, 150)
        ad = [random.randrange(N) for _ in range(L)]
        s = [random.randint(0, 1) for _ in range(N)]
        x = [random.randint(0, 1) for _ in range(L)]
        bad += frozen(s, ad, N, frozen(s, ad, N, x, 5), 6) != x
    print(f"  frozen addressing: W6 o W5 = id (arbitrary address sequences): {bad == 0}")
    fail = tot = 0
    for _ in range(2000):
        m = random.choice([1, 2, 3]); N = 1 << m; L = random.randint(4, 120)
        s = [random.randint(0, 1) for _ in range(N)]
        x = [random.randint(0, 1) for _ in range(L)]
        tot += 1
        fail += cxor(s, [0]*m, cxor(s, [0]*m, x, 5, m)[0], 6, m)[0] != x
    print(f"  UNFROZEN cxor^W6(s,pi,cxor^W5(s,pi,x)) != x: fails {fail}/{tot} "
          f"({100*fail/tot:.1f}%)")
    bad = tot = 0
    for _ in range(3000):
        m = random.choice([1, 2, 3]); N = 1 << m; L = random.randint(1, 120)
        pi = [random.randint(0, 1) for _ in range(m)]
        x = [random.randint(0, 1) for _ in range(L)]
        r5, addr = cxor([0]*N, pi, x, 5, m); r6, _ = cxor([0]*N, pi, x, 6, m)
        cnt = {}
        for j, a in enumerate(addr):
            v = cnt.get(a, 0); cnt[a] = v + 1
            tot += 1
            if r5[j] != r6[j] and v < 2:
                bad += 1
    print(f"  W5 and W6 agree at visit indices 0 and 1 of every class: {bad == 0} "
          f"({bad}/{tot} positions)")
    bad = 0
    for _ in range(2000):
        m = random.choice([1, 3]); N = 1 << m; L = 200
        s = [random.randint(0, 1) for _ in range(N)]
        s2 = [random.randint(0, 1) for _ in range(N)]
        x = [random.randint(0, 1) for _ in range(L)]
        u, ad = cxor(s, [0]*m, x, 5, m)
        bad += frozen(s2, ad, N, u, 6) != [x[j] ^ s[ad[j]] ^ s2[ad[j]] for j in range(L)]
    print(f"  seed mismatch gives x(k) + (s+s')[a_k] at every position: {bad == 0}")


def check_escape_hatches():
    """S11 -- the two ways to break affineness in the table."""
    def run_q(s, addr, N, x, upd, q):
        s = list(s); r = []
        for sym, a in zip(x, addr):
            c = s[a]; r.append((sym + c) % q); s[a] = upd(c, sym) % q
        return r

    def affine(upd, q, N=3, L=40, trials=300):
        for _ in range(trials):
            addr = [random.randrange(N) for _ in range(L)]
            x = [random.randrange(q) for _ in range(L)]
            s = [random.randrange(q) for _ in range(N)]
            t = [random.randrange(q) for _ in range(N)]
            lhs = run_q([(p + u) % q for p, u in zip(s, t)], addr, N, x, upd, q)
            z = run_q([0]*N, addr, N, x, upd, q)
            a = run_q(s, addr, N, x, upd, q); b = run_q(t, addr, N, x, upd, q)
            if lhs != [(p + u - v) % q for p, u, v in zip(a, b, z)]:
                return False
        return True
    print(f"  q=2, s <- y                      : affine in the table {affine(lambda c,y: y, 2)}")
    print(f"  q=4, s <- s*y   (scalar mult.)   : affine in the table {affine(lambda c,y: c*y, 4)}")
    print(f"  q=4, s <- s^2+y (nonlinear in s) : affine in the table {affine(lambda c,y: c*c+y, 4)}")

    def cxor_fb(s, pi, x, k, m):
        w = w_of(k); s = list(s); cxt = list(pi); r = []
        for sym in x:
            a = int("".join(map(str, cxt)), 2) if m else 0
            c = s[a]; o = sym ^ c; r.append(o); s[a] = w(c, sym)
            if m:
                cxt.pop(0); cxt.append(o)
        return r
    bad = set()
    for k in range(16):
        for _ in range(200):
            m = random.choice([1, 2]); N = 1 << m; L = random.randint(4, 40)
            x = [random.randint(0, 1) for _ in range(L)]
            s = [random.randint(0, 1) for _ in range(N)]
            t = [random.randint(0, 1) for _ in range(N)]
            lhs = cxor_fb([p ^ q for p, q in zip(s, t)], [0]*m, x, k, m)
            rhs = [a ^ b ^ c for a, b, c in zip(cxor_fb(s, [0]*m, x, k, m),
                                                cxor_fb(t, [0]*m, x, k, m),
                                                cxor_fb([0]*N, [0]*m, x, k, m))]
            if lhs != rhs:
                bad.add(k); break
    print(f"  address fed from the output: affineness fails for all 16 wirings: "
          f"{len(bad) == 16}")


def check_horizon():
    """S7/S15 -- the memory horizon has four regimes, indexed by (d1,d3)."""
    # Own generator, seeded to match gen_tables.emit_horizon, so that the numbers
    # here are the numbers printed in the paper's table regardless of run order.
    rng = random.Random(11)
    m, L = 4, 400; N = 1 << m; pi = [0] * m
    x = [rng.randint(0, 1) for _ in range(L)]

    def support(k, x):
        base, addr = cxor([0] * N, pi, x, k, m)
        sup = set()
        for a in range(N):
            e = [0] * N; e[a] = 1
            sup |= {j for j, (p, q) in enumerate(zip(cxor(e, pi, x, k, m)[0], base))
                    if p != q}
        return sup, addr

    # lambda_j = prod_{i<j} beta(y_i), beta = d1 + d3*y, position by position
    def predicted(k, x, addr):
        _, d1, _, d3 = anf(k)
        live = {a: 1 for a in range(N)}; pred = set()
        for j, a in enumerate(addr):
            if live[a]:
                pred.add(j)
            live[a] &= d1 ^ (d3 & x[j])
        return pred

    groups, bad = {}, []
    for k in range(16):
        sup, addr = support(k, x)
        if sup != predicted(k, x, addr):
            bad.append(k)
        groups.setdefault(anf(k)[1::2], []).append((k, len(sup)))
    print(f"  supp P matches lambda_j = prod beta(y_i) position by position, "
          f"all 16 wirings: {not bad}")
    NAME = {(0, 0): "first visit only", (1, 0): "never forgets",
            (0, 1): "survives a run of 1s", (1, 1): "survives a run of 0s"}
    for key in sorted(groups):
        ws = ", ".join(f"W{k}" for k, _ in groups[key])
        ns = sorted({n for _, n in groups[key]})
        print(f"    (d1,d3)={key}: {ws:<22} |supp P| = {ns[0] if len(ns)==1 else ns}"
              f"   {NAME[key]}")

    # The two intermediate regimes are tape-dependent; the two extremes are not.
    print("  tape bias moves the d3=1 regimes and not the d3=0 ones")
    print("  (the (0,0) row is |V|, which shrinks on a biased tape; its horizon "
          "does not move):")
    for bias in (0.1, 0.5, 0.9):
        xb = [1 if rng.random() < bias else 0 for _ in range(L)]
        row = {}
        for k in range(16):
            row.setdefault(anf(k)[1::2], set()).add(len(support(k, xb)[0]))
        nV = len(set(cxor([0] * N, pi, xb, 5, m)[1]))
        cells = "  ".join(f"{key}:{sorted(v)}" for key, v in sorted(row.items()))
        print(f"    P[y=1]={bias}: |V|={nV:2d}   {cells}")


def check_identifiability():
    """S10 -- how much tape it takes to tell two wirings apart."""
    # Own generator, seeded to match gen_tables.emit_identify.
    rng = random.Random(11)
    # (1) no wiring is visible at a first visit
    bad = tot = 0
    for m in (0, 1, 2, 3):
        N = 1 << m
        for L in range(1, 9):
            for t in product([0, 1], repeat=L):
                x = list(t)
                outs = [cxor([0] * N, [0] * m, x, k, m)[0] for k in range(16)]
                addr = cxor([0] * N, [0] * m, x, 5, m)[1]
                seen, firsts = set(), []
                for j, a in enumerate(addr):
                    if a not in seen:
                        seen.add(a); firsts.append(j)
                for k1 in range(16):
                    for k2 in range(k1 + 1, 16):
                        for j in firsts:
                            tot += 1
                            bad += outs[k1][j] != outs[k2][j]
    print(f"  no wiring is visible at a first visit: {bad == 0} "
          f"({bad}/{tot} disagreements, all 120 pairs, m<=3, tapes to length 8)")

    # (2)/(3) observational groups on a zero table: the partition is w(0,-),
    #         which the index convention puts in the top two bits, so it is k >> 2
    g = {}
    for k in range(16):
        g.setdefault((w_of(k)(0, 0), w_of(k)(0, 1)), []).append(k)
    ok = all(len({k >> 2 for k in v}) == 1 for v in g.values()) and len(g) == 4
    print(f"  the four observational groups are the index blocks k>>2: {ok}")
    for key in sorted(g):
        print(f"    w(0,0),w(0,1) = {key}: {['W%d' % k for k in g[key]]}")

    # (4) minimal separating tape for W5 vs W6, exhaustive
    print("  minimal tape separating W5 from W6 (zero table, pi = 0):")
    for m in range(0, 6):
        N = 1 << m; found = None
        for L in range(1, m + 6):
            for t in product([0, 1], repeat=L):
                x = list(t)
                if cxor([0]*N, [0]*m, x, 5, m)[0] != cxor([0]*N, [0]*m, x, 6, m)[0]:
                    found = (L, "".join(map(str, x))); break
            if found:
                break
        wit = [1] + [0] * (m + 2)
        sep = cxor([0]*N, [0]*m, wit, 5, m)[0] != cxor([0]*N, [0]*m, wit, 6, m)[0]
        print(f"    m={m}: length {found[0]} = m+3 ({found[0] == m + 3}), "
              f"first witness {found[1]}; 1 0^(m+2) separates: {sep}")

    # (5) median first divergence on random tapes -- a three-way birthday wait
    print("  first divergence of W5 and W6 on random tapes (200 trials):")
    for m in range(0, 9):
        N = 1 << m; L = max(64, 40 * N); pos = []
        for _ in range(200):
            x = [rng.randint(0, 1) for _ in range(L)]
            a = cxor([0]*N, [0]*m, x, 5, m)[0]; b = cxor([0]*N, [0]*m, x, 6, m)[0]
            d = [j for j, (p, q) in enumerate(zip(a, b)) if p != q]
            if d:
                pos.append(d[0])
        pos.sort()
        med = pos[len(pos) // 2]
        print(f"    m={m}: median {med:5d}   min {pos[0]:4d}   max {pos[-1]:5d}   "
              f"(2^(2m/3) = {2 ** (2 * m / 3):7.1f}, diverged {len(pos)}/200)")

    # a worked instance at m=4
    EX = "101001110111011111011111111111110"
    x = [int(c) for c in EX]
    r5 = "".join(map(str, cxor([0]*16, [0]*4, x, 5, 4)[0]))
    r6 = "".join(map(str, cxor([0]*16, [0]*4, x, 6, 4)[0]))
    d = [i for i, (p, q) in enumerate(zip(r5, r6)) if p != q]
    print(f"  worked instance, m=4, |x|={len(x)}: W5 and W6 first differ at "
          f"{d[0]}, Hamming distance {len(d)}")


def check_probe_recovery():
    """S14 -- the table is a linear unknown: one known (x,r) pair determines it."""
    def solve(pi, x, r, k, m):
        N = 1 << m
        beta, addr = cxor([0] * N, pi, x, k, m)
        cols = []
        for a in range(N):
            e = [0] * N; e[a] = 1
            cols.append([p ^ q for p, q in zip(cxor(e, pi, x, k, m)[0], beta)])
        rhs = [p ^ q for p, q in zip(r, beta)]
        aug = [[cols[a][i] for a in range(N)] + [rhs[i]] for i in range(len(rhs))]
        piv, rk = {}, 0
        for c in range(N):
            p = next((i for i in range(rk, len(aug)) if aug[i][c]), None)
            if p is None:
                continue
            aug[rk], aug[p] = aug[p], aug[rk]
            for i in range(len(aug)):
                if i != rk and aug[i][c]:
                    aug[i] = [u ^ v for u, v in zip(aug[i], aug[rk])]
            piv[c] = rk; rk += 1
        sol = [0] * N
        for c, i in piv.items():
            sol[c] = aug[i][N]
        return sol, set(addr), rk

    rng = random.Random(11)
    for k in (5, 6):
        ok = tot = 0
        for _ in range(300):
            m = rng.choice([1, 2, 3, 4]); N = 1 << m
            L = rng.randint(4, 200)
            s = [rng.randint(0, 1) for _ in range(N)]
            pi = [rng.randint(0, 1) for _ in range(m)]
            x = [rng.randint(0, 1) for _ in range(L)]
            r, _ = cxor(s, pi, x, k, m)
            sol, V, rk = solve(pi, x, r, k, m)
            tot += 1
            ok += all(sol[a] == s[a] for a in V) and rk == len(V)
        print(f"  W{k}: one known (x,r) pair recovers s on every visited address, "
              f"rank = |V|: {ok}/{tot}")
    print("  P and beta depend only on (pi,x), so no query to the machine is needed")


if __name__ == "__main__":
    random.seed(0)
    for fn in (check_index_and_lattice, check_causality_and_bijectivity,
               check_affine_in_table, check_support, check_not_affine_in_prefix,
               check_transfer, check_gray, check_duality_and_inverses,
               check_escape_hatches, check_horizon, check_identifiability,
               check_probe_recovery):
        print("=" * 76)
        print(f"[{fn.__name__}]  {fn.__doc__.splitlines()[0]}")
        fn()
