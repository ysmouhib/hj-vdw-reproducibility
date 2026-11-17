#!/usr/bin/env python3
"""Unrestricted one-hot van der Waerden SAT encoding with CryptoMiniSat."""
from __future__ import annotations

import argparse
import hashlib
import time

from cms_util import (
    add_at_least_one,
    add_first_occurrence,
    new_solver,
    one_hot_var,
    selected_colours,
    solve,
    write_json,
)


def progressions(N: int, k: int):
    for step in range(1, (N - 1) // (k - 1) + 1):
        for start in range(N - (k - 1) * step):
            yield tuple(start + j * step for j in range(k))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--symmetry", action="store_true")
    ap.add_argument("--time-limit", type=float)
    ap.add_argument("--confl-limit", type=int)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--output")
    args = ap.parse_args()
    r, k, N = args.r, args.k, args.N
    if r < 2 or k < 2 or N < k:
        ap.error("invalid r, k, or N")
    started = time.perf_counter()
    aps = list(progressions(N, k))
    solver = new_solver(args.time_limit, args.confl_limit, args.threads)
    add_at_least_one(solver, N, r)
    for ap_tuple in aps:
        for c in range(r):
            solver.add_clause([-one_hot_var(i, c, r) for i in ap_tuple])
    if args.symmetry:
        add_first_occurrence(solver, N, r)
    print(
        f"BUILT N={N} r={r} k={k} variables={N*r} APs={len(aps)} "
        f"clauses={N+r*len(aps)} fixed_positions=0 symmetry={args.symmetry}",
        flush=True,
    )
    status, model = solve(solver)
    print(f"RESULT {status} seconds={time.perf_counter()-started:.3f}")
    if status != "SAT":
        return 1 if status == "UNKNOWN" else 0
    colours = selected_colours(model, N, r)
    bad = sum(len({colours[i] for i in ap_tuple}) == 1 for ap_tuple in aps)
    if bad:
        raise RuntimeError(f"model failed direct check: {bad} monochromatic APs")
    print(f"VERIFIED APs={len(aps)} bad=0")
    if args.output:
        colouring = "".join(chr(ord("a") + c) for c in colours)
        write_json(args.output, {
            "kind": "van-der-Waerden-coloring",
            "r": r,
            "k": k,
            "N": N,
            "positions": "1,...,N",
            "alphabet": "a=0, b=1, ...",
            "coloring": colouring,
            "sha256_coloring": hashlib.sha256(colouring.encode()).hexdigest(),
            "monochromatic_progressions": 0,
        })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
