#!/usr/bin/env python3
"""Original-style global vdW CDCL run: warm phases, no fixed positions."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from vdw_cms import progressions


def read_warm(path: str) -> tuple[int, int, int, list[int]]:
    data = json.loads(Path(path).read_text())
    text = data.get("coloring", data.get("palette"))
    if text is None:
        raise ValueError("warm file needs coloring or palette")
    if all("a" <= c <= "z" for c in text):
        colours = [ord(c) - ord("a") for c in text]
    else:
        colours = [int(c, 36) for c in text]
    return data["r"], data["k"], data["N"], colours


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--warm", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--solver", default="cadical300")
    ap.add_argument("--budget", type=int, default=1_000_000, help="conflict budget")
    ap.add_argument("--symmetry", action="store_true")
    args = ap.parse_args()
    try:
        from pysat.solvers import Solver
    except ImportError as exc:
        raise SystemExit("Install python-sat for this optional CDCL driver.") from exc
    r, k, N, warm = read_warm(args.warm)
    if len(warm) != N or any(not 0 <= c < r for c in warm):
        raise ValueError("invalid warm colouring")
    started = time.perf_counter()
    var = lambda i, c: i * r + c + 1
    aps = list(progressions(N, k))
    with Solver(name=args.solver) as solver:
        for i in range(N):
            solver.add_clause([var(i, c) for c in range(r)])
        if args.symmetry:
            for i in range(N):
                for c in range(1, r):
                    solver.add_clause([-var(i, c)] + [var(h, c - 1) for h in range(i)])
        for ap_tuple in aps:
            for c in range(r):
                solver.add_clause([-var(i, c) for i in ap_tuple])
        # These are phase preferences only: no warm value is fixed by a clause.
        try:
            solver.set_phases([var(i, c) if warm[i] == c else -var(i, c)
                               for i in range(N) for c in range(r)])
        except (NotImplementedError, AttributeError):
            print("PHASES unsupported by selected solver; continuing without them")
        print(f"BUILT vars={N*r} APs={len(aps)} symmetry={args.symmetry}", flush=True)
        solver.conf_budget(args.budget)
        result = solver.solve_limited()
        status = {True: "SAT", False: "UNSAT", None: "UNKNOWN"}[result]
        print(f"RESULT {status} conflict_budget={args.budget} seconds={time.perf_counter()-started:.3f}")
        if result is not True:
            return 1 if result is None else 0
        model = set(solver.get_model())
    colours = [next(c for c in range(r) if var(i, c) in model) for i in range(N)]
    bad = sum(len({colours[i] for i in ap_tuple}) == 1 for ap_tuple in aps)
    if bad:
        raise RuntimeError(f"model failed direct check: {bad} monochromatic APs")
    colouring = "".join(chr(ord("a") + c) for c in colours)
    Path(args.output).write_text(json.dumps({
        "kind": "van-der-Waerden-coloring", "r": r, "k": k, "N": N,
        "positions": "1,...,N", "alphabet": "a=0, b=1, ...",
        "coloring": colouring,
        "sha256_coloring": hashlib.sha256(colouring.encode()).hexdigest(),
        "monochromatic_progressions": 0,
    }, indent=2) + "\n")
    print(f"VERIFIED APs={len(aps)} bad=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
