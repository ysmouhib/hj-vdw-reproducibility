"""Small, explicit helpers for the CryptoMiniSat 5.18.0 Python interface."""
from __future__ import annotations

import json
from pathlib import Path


def new_solver(time_limit=None, confl_limit=None, threads=1, verbose=0):
    try:
        from pycryptosat import Solver
    except ImportError as exc:
        raise SystemExit(
            "pycryptosat is required for search, not verification: "
            "python3 -m pip install pycryptosat==5.18.0"
        ) from exc
    kwargs = {"threads": threads, "verbose": verbose}
    if time_limit is not None:
        kwargs["time_limit"] = time_limit
    if confl_limit is not None:
        kwargs["confl_limit"] = confl_limit
    return Solver(**kwargs)


def solve(solver, assumptions=None):
    """Return SAT/UNSAT/UNKNOWN; a budget-limited None is never called UNSAT."""
    sat, model = solver.solve() if assumptions is None else solver.solve(assumptions=assumptions)
    if sat is None:
        return "UNKNOWN", None
    return ("SAT", model) if sat else ("UNSAT", None)


def one_hot_var(i: int, colour: int, r: int) -> int:
    return i * r + colour + 1


def add_at_least_one(solver, cells: int, r: int) -> None:
    for i in range(cells):
        solver.add_clause([one_hot_var(i, c, r) for c in range(r)])


def add_at_most_one(solver, cells: int, r: int) -> None:
    for i in range(cells):
        for c in range(r):
            for d in range(c + 1, r):
                solver.add_clause([-one_hot_var(i, c, r), -one_hot_var(i, d, r)])


def add_first_occurrence(solver, cells: int, r: int) -> None:
    """Break colour renaming: colour c requires an earlier occurrence of c-1."""
    for i in range(cells):
        for c in range(1, r):
            solver.add_clause(
                [-one_hot_var(i, c, r)]
                + [one_hot_var(h, c - 1, r) for h in range(i)]
            )


def selected_colours(model, cells: int, r: int) -> list[int]:
    """Extract an ordinary colouring; multiple true colours are harmless here."""
    out = []
    for i in range(cells):
        chosen = [c for c in range(r) if model[one_hot_var(i, c, r)]]
        if not chosen:
            raise ValueError(f"cell {i} has no selected colour")
        out.append(chosen[0])
    return out


def add_sequential_at_most(solver, literals: list[int], bound: int, next_var: int):
    """Equational sequential counter for sum(literals) <= bound.

    s[i,j] means that at least j+1 inputs among literals[0:i+1] are true.
    The final flag s[n-1,bound] (at least bound+1) is forced false.
    Returns the next unused DIMACS variable number.
    """
    n = len(literals)
    if bound < 0:
        raise ValueError("negative cardinality bound")
    if bound >= n:
        return next_var
    if bound == 0:
        for lit in literals:
            solver.add_clause([-lit])
        return next_var
    # j=0,...,bound; the last column means at least bound+1.
    s = [[next_var + i * (bound + 1) + j for j in range(bound + 1)]
         for i in range(n)]
    next_var += n * (bound + 1)
    solver.add_clause([-s[0][0], literals[0]])
    solver.add_clause([-literals[0], s[0][0]])
    for j in range(1, bound + 1):
        solver.add_clause([-s[0][j]])
    for i in range(1, n):
        x = literals[i]
        # j=0: s[i,0] <-> s[i-1,0] OR x.
        solver.add_clause([-s[i][0], s[i - 1][0], x])
        solver.add_clause([-s[i - 1][0], s[i][0]])
        solver.add_clause([-x, s[i][0]])
        for j in range(1, bound + 1):
            prev = s[i - 1][j]
            prev_lower = s[i - 1][j - 1]
            cur = s[i][j]
            # cur <-> prev OR (x AND prev_lower).
            solver.add_clause([-cur, prev, x])
            solver.add_clause([-cur, prev, prev_lower])
            solver.add_clause([-prev, cur])
            solver.add_clause([-x, -prev_lower, cur])
    solver.add_clause([-s[n - 1][bound]])
    return next_var


def write_json(path: str | Path, data) -> None:
    Path(path).write_text(json.dumps(data, indent=2) + "\n")
