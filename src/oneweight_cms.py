#!/usr/bin/env python3
"""Exact CryptoMiniSat search for a one-weight Hales--Jewett palette."""
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
from simplex_cms import compositions


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", required=True, help="comma-separated integer weights")
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--modulus", type=int, help="optionally identify levels modulo m")
    ap.add_argument("--symmetry", action="store_true")
    ap.add_argument("--time-limit", type=float)
    ap.add_argument("--confl-limit", type=int)
    ap.add_argument("--output")
    args = ap.parse_args()
    omega = tuple(int(x) for x in args.omega.split(","))
    t, n, r = len(omega), args.n, args.r
    if t < 2 or n < 1 or r < 2:
        ap.error("invalid weights, n, or r")

    level_ids = {}
    levels = []

    def level_id(value):
        if args.modulus:
            value %= args.modulus
        if value not in level_ids:
            level_ids[value] = len(levels)
            levels.append(value)
        return level_ids[value]

    for u in compositions(n, t):
        level_id(sum(a * b for a, b in zip(omega, u)))
    homothets = set()
    for k in range(1, n + 1):
        for u in compositions(n - k, t):
            b = sum(a * b for a, b in zip(omega, u))
            homothets.add(tuple(level_id(b + k * omega[a]) for a in range(t)))
    homothets = sorted(homothets)
    started = time.perf_counter()
    solver = new_solver(args.time_limit, args.confl_limit)

    if r == 2:
        # One Boolean per level: false/true are the two colours.
        for homothet in homothets:
            solver.add_clause([i + 1 for i in homothet])
            solver.add_clause([-(i + 1) for i in homothet])
        # Harmless global complement symmetry.
        if 0 in level_ids:
            solver.add_clause([-(level_ids[0] + 1)])
    else:
        add_at_least_one(solver, len(levels), r)
        for homothet in homothets:
            for c in range(r):
                solver.add_clause([-one_hot_var(i, c, r) for i in homothet])
        if args.symmetry:
            add_first_occurrence(solver, len(levels), r)

    print(f"BUILT levels={len(levels)} homothets={len(homothets)} r={r}", flush=True)
    status, model = solve(solver)
    print(f"RESULT {status} seconds={time.perf_counter()-started:.3f}")
    if status != "SAT":
        return 1 if status == "UNKNOWN" else 0
    if r == 2:
        colours = [int(bool(model[i + 1])) for i in range(len(levels))]
    else:
        colours = selected_colours(model, len(levels), r)
    bad = sum(len({colours[i] for i in h}) == 1 for h in homothets)
    if bad:
        raise RuntimeError(f"model failed direct check: {bad} monochromatic homothets")
    print(f"VERIFIED levels={len(levels)} homothets={len(homothets)} bad=0")
    if args.output:
        alphabet = "0123456789abcdefghijklmnopqrstuvwxyz"
        palette = "".join(alphabet[c] for c in colours)
        write_json(args.output, {
            "kind": "one-weight-coloring",
            "omega": omega,
            "r": r,
            "n": n,
            "modulus": args.modulus,
            "level_order": levels,
            "palette": palette,
            "sha256_palette": hashlib.sha256(palette.encode()).hexdigest(),
            "monochromatic_homothets": 0,
        })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
