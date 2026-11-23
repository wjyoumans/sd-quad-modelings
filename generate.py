#!/usr/bin/env python3
"""Planted syndrome decoding instances at the Gilbert-Varshamov weight.

Usage: generate.py n k seed

Writes instances/SD_<n>_<k>_<seed>, next to this script:
- H: a uniformly random (n-k) x n matrix over F_2, resampled until it has
  full rank (not in systematic form);
- t = t_GV(n, k), the largest t <= n/2 with binom(n, t) <= 2^(n-k);
- e: a uniformly random vector of weight t, and s = H e.

The file has the sections "# n", "# k", "# seed", "# w" (= t),
"# H (each line is a row of H)", "# s^transpose" and "# solution" (= e).
The output depends only on (n, k, seed): the random generator is seeded with
the string "SD_GV n k seed", and it only uses getrandbits and sample, whose
output for a given seed is stable across Python 3 versions.
"""

import os
import random
import sys
from math import comb


def t_gv(n, k):
    return max(t for t in range(n // 2 + 1) if comb(n, t) <= 2 ** (n - k))


# The values the experiments rely on.
assert all(t_gv(n, n // 2) == t for n, t in {16: 2, 20: 2, 24: 3, 28: 3, 32: 4, 36: 4, 40: 5,
                                            44: 5, 48: 6, 56: 7, 60: 7}.items())
assert all(t_gv(n, n // 4) == n // 4 for n in range(16, 45, 4))


def rank(rows):
    """Rank over F_2 of the rows, given as integers."""
    pivots = {}           # leading bit -> a row with that leading bit
    for r in rows:
        while r:
            lead = r.bit_length() - 1
            if lead not in pivots:
                pivots[lead] = r
                break
            r ^= pivots[lead]
    return len(pivots)


def generate(n, k, seed):
    rng = random.Random(f"SD_GV {n} {k} {seed}")
    m, t = n - k, t_gv(n, k)
    while True:
        H = [[rng.getrandbits(1) for _ in range(n)] for _ in range(m)]
        if rank([int("".join(map(str, row)), 2) for row in H]) == m:
            break
    support = set(rng.sample(range(n), t))
    e = [1 if j in support else 0 for j in range(n)]
    s = [sum(h * x for h, x in zip(row, e)) % 2 for row in H]
    return H, s, e, t


def bits(v):
    return "".join(map(str, v))


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    n, k, seed = map(int, sys.argv[1:])
    if not 1 <= k < n:
        sys.exit("need 1 <= k < n")
    H, s, e, t = generate(n, k, seed)
    text = (f"# n\n{n}\n# k\n{k}\n# seed\n{seed}\n# w\n{t}\n"
            "# H (each line is a row of H)\n" + "".join(bits(row) + "\n" for row in H)
            + f"# s^transpose\n{bits(s)}\n# solution\n{bits(e)}\n")
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "instances")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, f"SD_{n}_{k}_{seed}"), "w") as f:
        f.write(text)


if __name__ == "__main__":
    main()
