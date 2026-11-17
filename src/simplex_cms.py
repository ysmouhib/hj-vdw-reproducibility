#!/usr/bin/env python3
"""CryptoMiniSat encoding for symmetric/bracket Hales--Jewett colourings."""
from __future__ import annotations

import argparse
import hashlib
import time
from pathlib import Path

from cms_util import (
    add_at_least_one,
    add_at_most_one,
    add_first_occurrence,
    add_sequential_at_most,
    new_solver,
    one_hot_var,
    selected_colours,
    solve,
    write_json,
)


def compositions(total: int, parts: int, prefix=()):
    if parts == 1:
        yield prefix + (total,)
    else:
        for x in range(total + 1):
            yield from compositions(total - x, parts - 1, prefix + (x,))


def make_instance(t: int, n: int, max_k: int):
    cells = list(compositions(n, t))
    index = {v: i for i, v in enumerate(cells)}
    corners = []
    for k in range(1, min(n, max_k) + 1):
        for u in compositions(n - k, t):
            corner = []
            for a in range(t):
                v = list(u)
                v[a] += k
                corner.append(index[tuple(v)])
            corners.append(tuple(corner))
    return cells, corners


def bad_corners(colours, corners) -> int:
    return sum(len({colours[i] for i in corner}) == 1 for corner in corners)


def read_warm(path: str, cells: list[tuple[int]]) -> list[int]:
    import json

    data = json.loads(Path(path).read_text())
    if "rows" in data:
        if len(cells[0]) != 3:
            raise ValueError("row input is only for t=3")
        rows = data["rows"]
        return [int(rows[v[0]][v[1]]) for v in cells]
    if "colouring" in data:
        values = [int(x, 36) for x in data["colouring"]]
    elif "colors" in data:
        values = data["colors"]
    else:
        raise ValueError("warm JSON needs rows, colouring, or colors")
    if len(values) != len(cells):
        raise ValueError("warm colouring has the wrong number of cells")
    return values


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--t", type=int, required=True)
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--max-k", type=int, help="restrict lines to at most this many active coordinates")
    ap.add_argument("--symmetry", action="store_true", help="complete first-occurrence colour breaking")
    ap.add_argument("--fix-first", action="store_true", help="fix cell 0 to colour 0 (original T21 calibration)")
    ap.add_argument("--warm", help="warm JSON, required with --radius")
    ap.add_argument("--radius", type=int, help="exact Hamming neighbourhood around --warm")
    ap.add_argument("--time-limit", type=float)
    ap.add_argument("--confl-limit", type=int)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--verbose", type=int, default=0)
    ap.add_argument("--output")
    args = ap.parse_args()
    if args.t < 2 or args.r < 2 or args.n < 1:
        ap.error("require t,r >= 2 and n >= 1")
    if args.radius is not None and not args.warm:
        ap.error("--radius requires --warm")
    if args.radius is not None and args.symmetry:
        ap.error("do not combine colour symmetry with a warm-centred Hamming ball")

    started = time.perf_counter()
    cells, corners = make_instance(args.t, args.n, args.max_k or args.n)
    r = args.r
    solver = new_solver(args.time_limit, args.confl_limit, args.threads, args.verbose)
    add_at_least_one(solver, len(cells), r)
    for corner in corners:
        for c in range(r):
            solver.add_clause([-one_hot_var(i, c, r) for i in corner])
    if args.symmetry:
        add_first_occurrence(solver, len(cells), r)
    if args.fix_first:
        solver.add_clause([one_hot_var(0, 0, r)])
    warm = None
    if args.radius is not None:
        warm = read_warm(args.warm, cells)
        if any(not 0 <= c < r for c in warm):
            ap.error("warm colour outside palette")
        # Exactly one colour per cell makes Hamming distance unambiguous.
        add_at_most_one(solver, len(cells), r)
        changes = [-one_hot_var(i, warm[i], r) for i in range(len(cells))]
        add_sequential_at_most(solver, changes, args.radius, len(cells) * r + 1)

    print(
        f"BUILT cells={len(cells)} corners={len(corners)} variables_at_least={len(cells)*r} "
        f"symmetry={args.symmetry} radius={args.radius}",
        flush=True,
    )
    status, model = solve(solver)
    print(f"RESULT {status} seconds={time.perf_counter()-started:.3f}", flush=True)
    if status == "UNSAT" and args.radius is not None:
        print("SCOPE: UNSAT only inside the stated Hamming neighbourhood; not a global bound.")
    if status != "SAT":
        return 1 if status == "UNKNOWN" else 0
    colours = selected_colours(model, len(cells), r)
    bad = bad_corners(colours, corners)
    if bad:
        raise RuntimeError(f"solver model failed direct check: {bad} bad corners")
    print(f"VERIFIED corners={len(corners)} bad=0")
    if args.output:
        data = {
            "kind": "symmetric-Hales-Jewett-simplex-certificate",
            "t": args.t,
            "r": r,
            "n": args.n,
            "max_k": args.max_k or args.n,
            "monochromatic_corners": 0,
            "corner_tuples_checked": len(corners),
            "cell_order": "recursive lexicographic weak compositions",
            "colors": colours,
        }
        if args.t == 3:
            data["rows"] = [
                "".join(str(colours[cells.index((a, b, args.n-a-b))])
                        for b in range(args.n-a+1))
                for a in range(args.n+1)
            ]
            data["sha256_rows"] = hashlib.sha256("".join(data["rows"]).encode()).hexdigest()
        write_json(args.output, data)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
