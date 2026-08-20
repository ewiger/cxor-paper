"""
CTA update-rule study: W5 (code, A.1) vs W6 (paper p.8 step 10 / p.12 eq. 1).

  W5  "last value"     S[a] <- y(k)
  W6  "running parity" S[a] <- S[a] ^ y(k)   (== the emitted bit r(k))

Both: r(k) = y(k) ^ S[a_k],  a_k = y(k-m .. k-1)  (window EXCLUDES current bit;
that ordering, not the paper's steps 6-7, is what reproduces the p.9 example).
"""
from itertools import product
import random

EXAMPLE_IN  = "101001110111011111011111111111110"
EXAMPLE_OUT = "101001110111000011100001000000001"   # as printed on p.9 / in the A.1 docstring


def fwd(y, m, x=None, rule="W5", order="exclude"):
    """Encode. Returns (output bits, address sequence).

    order="exclude": address is y(k-m..k-1), read BEFORE the context shift  (code A.1)
    order="include": address is y(k-m+1..k), read AFTER  the context shift  (paper steps 6-7)
    """
    n = 1 << m
    S = list(x) if x else [0] * n
    cxt, out, addr = [0] * m, [], []
    for b in y:
        if order == "include" and m:
            cxt.pop(0); cxt.append(b)
        a = int("".join(map(str, cxt)), 2) if m else 0
        addr.append(a)
        s = S[a]; r = b ^ s
        out.append(r)
        S[a] = b if rule == "W5" else r
        if order == "exclude" and m:
            cxt.pop(0); cxt.append(b)
    return out, addr


def inv(r, m, x=None, rule="W5"):
    """Decode (addressing rebuilt from recovered bits). Inverse of fwd for both rules."""
    n = 1 << m
    S = list(x) if x else [0] * n
    cxt, out = [0] * m, []
    for rk in r:
        a = int("".join(map(str, cxt)), 2) if m else 0
        b = rk ^ S[a]
        out.append(b)
        S[a] = b if rule == "W5" else rk
        if m:
            cxt.pop(0); cxt.append(b)
    return out


def frozen(tape, addr, naddr, x=None, rule="W5"):
    """Run with an EXTERNALLY supplied address sequence (no window, no y-dependence)."""
    S = list(x) if x else [0] * naddr
    out = []
    for b, a in zip(tape, addr):
        s = S[a]; r = b ^ s
        out.append(r)
        S[a] = b if rule == "W5" else r
    return out


def visited(y, m):
    cxt, s = [0] * m, set()
    for b in y:
        s.add(int("".join(map(str, cxt)), 2) if m else 0)
        if m:
            cxt.pop(0); cxt.append(b)
    return s


def linear_part(y, m, rule):
    """Return (beta, columns) with cxor(x,y) = beta ^ sum_a x[a]*col_a."""
    n = 1 << m
    base, _ = fwd(y, m, [0] * n, rule)
    cols = []
    for a in range(n):
        e = [0] * n; e[a] = 1
        cols.append([p ^ q for p, q in zip(fwd(y, m, e, rule)[0], base)])
    return base, cols


# ----------------------------------------------------------------------------- checks

def check_example():
    y = [int(c) for c in EXAMPLE_IN]
    r1 = "".join(map(str, fwd(y, 4, rule="W5")[0]))
    r2 = "".join(map(str, fwd(y, 4, rule="W6")[0]))
    d = [i for i, (p, q) in enumerate(zip(r1, r2)) if p != q]
    print("p.9 example, m=4")
    print("  W5        ", r1, "MATCH" if r1 == EXAMPLE_OUT else "differs")
    print("  W6        ", r2, "MATCH" if r2 == EXAMPLE_OUT else "differs")
    print(f"  first divergence at index {d[0]}, Hamming distance {len(d)}/{len(y)}")


def check_shortest_divergence(ms=(1, 2, 4), maxn=14):
    print("shortest y (x=0) with W5(0,y) != W6(0,y)")
    for m in ms:
        for n in range(1, maxn):
            hit = next((y for y in product([0, 1], repeat=n)
                        if fwd(list(y), m, rule="W5")[0] != fwd(list(y), m, rule="W6")[0]), None)
            if hit:
                print(f"  m={m}: len={n}  y={''.join(map(str,hit))}")
                break


def check_bijective(ms=(1, 2, 3), maxn=10, trials=30):
    print("bijectivity in y (exhaustive over {0,1}^n, random x)")
    for m in ms:
        for rule in ("W5", "W6"):
            ok = True
            for _ in range(trials):
                x = [random.randint(0, 1) for _ in range(1 << m)]
                for n in range(1, maxn + 1):
                    if len({tuple(fwd(list(y), m, x, rule)[0])
                            for y in product([0, 1], repeat=n)}) != 2 ** n:
                        ok = False
            print(f"  m={m} {rule}: bijective -> {ok}")
    print("round-trip inv(fwd(y)) == y")
    for rule in ("W5", "W6"):
        ok = all(inv(fwd(y, m, x, rule)[0], m, x, rule) == y
                 for m, x, y in ((m, [random.randint(0, 1) for _ in range(1 << m)],
                                  [random.randint(0, 1) for _ in range(random.randint(1, 200))])
                                 for m in random.choices([1, 2, 4], k=1000)))
        print(f"  {rule}: {ok}")


def check_frozen_inverse():
    """W6_frozen(x, W5(x,y)) == y  -- uniform in m, and in fact in the addressing."""
    print("frozen-addressing inverse, W6 o W5 = id")
    for m in (1, 2, 3):
        for n in (1, 5, 10):
            nadr = 1 << m; bad = tot = 0
            for x in product([0, 1], repeat=nadr):
                for y in product([0, 1], repeat=n):
                    y = list(y); tot += 1
                    u, ad = fwd(y, m, list(x), "W5")
                    bad += frozen(u, ad, nadr, list(x), "W6") != y
            print(f"  m={m} n={n}: {tot} (x,y) pairs, failures={bad}")
    for m in (4, 6, 8, 10, 12):
        nadr = 1 << m; bad = 0
        for _ in range(200):
            x = [random.randint(0, 1) for _ in range(nadr)]
            y = [random.randint(0, 1) for _ in range(random.randint(1, 5000))]
            u, ad = fwd(y, m, x, "W5")
            bad += frozen(u, ad, nadr, x, "W6") != y
        print(f"  m={m}: 200 random (x,y), |y|<=5000, failures={bad}")
    bad = 0
    for _ in range(2000):
        nadr = random.randint(1, 20); L = random.randint(1, 300)
        ad = [random.randrange(nadr) for _ in range(L)]
        x = [random.randint(0, 1) for _ in range(nadr)]
        y = [random.randint(0, 1) for _ in range(L)]
        bad += frozen(frozen(y, ad, nadr, x, "W5"), ad, nadr, x, "W6") != y
    print(f"  arbitrary address sequences (not from any y): failures={bad}/2000")


def check_seed_mismatch():
    """Decoding with x2 != x offsets each class by (x^x2)[a] at EVERY position."""
    bad = 0
    for _ in range(2000):
        m = random.choice([1, 3, 5]); nadr = 1 << m
        x = [random.randint(0, 1) for _ in range(nadr)]
        x2 = [random.randint(0, 1) for _ in range(nadr)]
        y = [random.randint(0, 1) for _ in range(300)]
        u, ad = fwd(y, m, x, "W5")
        v = frozen(u, ad, nadr, x2, "W6")
        bad += v != [y[k] ^ x[ad[k]] ^ x2[ad[k]] for k in range(len(y))]
    print(f"seed-mismatch error == (x^x2)[a_k] at every k: {bad == 0} ({bad}/2000 deviations)")


def check_support_rank():
    """rank(P_y) <= 2^m under both; support finite under W5, cofinite under W6."""
    y = [random.randint(0, 1) for _ in range(300)]
    print("linear part of cxor_m(.,y), m=4, |y|=300")
    for rule in ("W5", "W6"):
        _, cols = linear_part(y, 4, rule)
        nz = [a for a, c in enumerate(cols) if any(c)]
        sup = {k for c in cols for k, v in enumerate(c) if v}
        print(f"  {rule}: rank={len(nz)}  |supp P_y|={len(sup)} of {len(y)}")


def check_ordering():
    """The window convention is forced: paper steps 6-7 (include) are NOT injective.

    Also resolves the update rule -- exactly one of the four readings reproduces p.9.
    """
    y = [int(c) for c in EXAMPLE_IN]
    print("p.9 example, m=4 -- all four readings of the spec")
    for order in ("exclude", "include"):
        for rule in ("W5", "W6"):
            got = "".join(map(str, fwd(y, 4, rule=rule, order=order)[0]))
            tag = "MATCH" if got == EXAMPLE_OUT else "differs"
            print(f"  order={order:8s} rule={rule}  {tag:8s} {got}")
    print("injectivity of y -> CTA(0,y), exhaustive m<=4, n<=12")
    for order in ("exclude", "include"):
        for rule in ("W5", "W6"):
            bad = None
            for m in (1, 2, 3, 4):
                for n in range(1, 13):
                    seen = {}
                    for t in product([0, 1], repeat=n):
                        k = tuple(fwd(list(t), m, rule=rule, order=order)[0])
                        if k in seen and bad is None:
                            bad = (m, "".join(map(str, seen[k])), "".join(map(str, t)))
                        seen[k] = t
            if bad:
                print(f"  order={order:8s} rule={rule}: NOT injective -- "
                      f"m={bad[0]}: y={bad[1]} and y'={bad[2]} collide")
            else:
                print(f"  order={order:8s} rule={rule}: injective")


def check_guard_and_table():
    """`if diff:` in A.1 is a no-op; `2 << m` double-allocates the table."""
    def cta(y, m, guard):
        S = [0] * (2 << m); cxt, out, used = [0] * m, [], set()
        for b in y:
            i = int("".join(map(str, cxt)), 2); used.add(i)
            s = S[i]; d = s ^ b; out.append(d)
            if not guard or d:
                S[i] = b
            cxt.pop(0); cxt.append(b)
        return out, used
    ys = [[random.randint(0, 1) for _ in range(200)] for _ in range(500)]
    print("`if diff:` guard is a no-op:",
          all(cta(y, 4, True)[0] == cta(y, 4, False)[0] for y in ys))
    used = set().union(*(cta(y, 4, True)[1] for y in ys))
    print(f"table `2 << 4` = {2<<4} cells; addresses ever used = {len(used)} (max {max(used)})")


if __name__ == "__main__":
    random.seed(0)
    for f in (check_example, check_ordering, check_shortest_divergence, check_bijective,
              check_frozen_inverse, check_seed_mismatch, check_support_rank,
              check_guard_and_table):
        print("=" * 72); f()