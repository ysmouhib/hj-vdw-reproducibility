#!/usr/bin/env python3
"""Independently verify every certificate in this repository (stdlib only)."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("ascii")).hexdigest()


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def prime_factors(n: int) -> list[int]:
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        out.append(n)
    return out


def is_primitive_root(p: int, g: int) -> bool:
    return 1 <= g < p and all(
        pow(g, (p - 1) // q, p) != 1 for q in prime_factors(p - 1)
    )


def power_character(p: int, r: int, g: int) -> list[int | None]:
    """chi[x] = ind_g(x) mod r, with chi[0] = None."""
    chi: list[int | None] = [None] * p
    x = 1
    for exponent in range(p - 1):
        if chi[x] is not None:
            raise ValueError("g does not enumerate all nonzero residues")
        chi[x] = exponent % r
        x = x * g % p
    if x != 1 or any(v is None for v in chi[1:]):
        raise ValueError("invalid primitive-root orbit")
    return chi


def verify_algebraic(path: Path) -> tuple[int, int, int]:
    doc = json.loads(path.read_text())
    records = doc["records"]
    total_translates = 0
    table_records = 0
    for rec in records:
        t, r, p, g = (rec[k] for k in ("t", "r", "p", "g"))
        S = rec["S"]
        assert len(S) == t and len(set(S)) == t
        assert all(type(s) is int and 0 <= s < p for s in S)
        assert is_prime(p), f"{p} is not prime"
        assert (p - 1) % r == 0
        assert is_primitive_root(p, g), f"g={g} is not primitive modulo {p}"
        assert len({s % 2 for s in S}) == 2, "mixed parity is required for p+1"
        chi = power_character(p, r, g)
        bad = []
        for y in range(p):
            colours = {chi[(y + s) % p] for s in S if (y + s) % p != 0}
            if len(colours) < 2:
                bad.append(y)
        assert not bad, f"(t,r,p)=({t},{r},{p}): failing translates {bad[:5]}"
        # The exact period-2p palette used by the article, including 0 and p.
        palette_chars = []
        for x in range(2 * p):
            residue = x % p
            if residue:
                colour = chi[residue]
            else:
                colour = 0 if x == 0 else 1
            palette_chars.append(str(colour))
        palette_hash = sha256_text("".join(palette_chars))
        if "palette_sha256" in rec:
            assert palette_hash == rec["palette_sha256"], f"palette hash: {p}"
        assert rec["HJ"] == p + 1 and rec["kappa_ow"] == p - 1
        assert rec.get("R_bad", 0) == 0
        total_translates += p
        table_records += rec.get("status") == "article_table"
        print(f"OK algebraic (t,r)=({t},{r}) p={p}: {p} translates, palette hash OK")
    return len(records), table_records, total_translates


def verify_simplex(path: Path) -> tuple[int, int]:
    rec = json.loads(path.read_text())
    t, r, n, rows = (rec[k] for k in ("t", "r", "n", "rows"))
    assert t == 3, "this compact row format is for t=3"
    assert len(rows) == n + 1
    for a, row in enumerate(rows):
        assert len(row) == n - a + 1
        assert all("0" <= c < str(r) for c in row)
    joined = "".join(rows)
    if "sha256_rows" in rec:
        assert sha256_text(joined) == rec["sha256_rows"], f"row hash: {path.name}"
    bad = 0
    corners = 0
    # Coordinates are (a,b,n-a-b); a line's three types differ by k*e_i.
    for k in range(1, n + 1):
        for a in range(n - k + 1):
            for b in range(n - k - a + 1):
                colours = (rows[a + k][b], rows[a][b + k], rows[a][b])
                corners += 1
                bad += colours[0] == colours[1] == colours[2]
    assert bad == 0, f"{path.name}: {bad} monochromatic corners"
    expected_cells = (n + 1) * (n + 2) // 2
    assert len(joined) == expected_cells
    if "cells" in rec:
        assert rec["cells"] == expected_cells
    if "corner_triples_checked" in rec:
        assert rec["corner_triples_checked"] == corners
    if "monochromatic_corners" in rec:
        assert rec["monochromatic_corners"] == 0
    print(f"OK simplex HJ(3,{r})>={n + 1}: {expected_cells} cells, {corners} corners")
    return expected_cells, corners


def verify_vdw(path: Path) -> int:
    rec = json.loads(path.read_text())
    r, k, N, s = (rec[k] for k in ("r", "k", "N", "coloring"))
    assert k == 3 and len(s) == N
    colours = [ord(c) - ord("a") for c in s]
    assert all(0 <= c < r for c in colours)
    if "sha256_coloring" in rec:
        assert sha256_text(s) == rec["sha256_coloring"], f"colouring hash: {path.name}"
    progressions = 0
    for start in range(N - 2):
        for step in range(1, (N - 1 - start) // 2 + 1):
            progressions += 1
            if colours[start] == colours[start + step] == colours[start + 2 * step]:
                raise AssertionError(
                    f"{path.name}: monochromatic AP at 1-based positions "
                    f"{start + 1},{start + step + 1},{start + 2 * step + 1}"
                )
    print(f"OK vdW W({r},3)>{N}: {progressions} progressions")
    return progressions


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    repository_root = script_dir.parent
    # Also support an arXiv ancillary layout: anc/verify_all.py + anc/certificates/.
    default_root = repository_root if (repository_root / "certificates").is_dir() else script_dir
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    args = parser.parse_args()
    root = args.root.resolve()
    cert = root / "certificates"
    manifest_path = cert / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for entry in manifest["files"]:
        path = cert / entry["path"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == entry["sha256"], f"file integrity: {entry['path']}"
    print(f"OK manifest: {len(manifest['files'])} certificate files")

    algebraic = verify_algebraic(cert / "algebraic.json")
    simplex_totals = [0, 0]
    for path in sorted((cert / "simplex").glob("*.json")):
        cells, corners = verify_simplex(path)
        simplex_totals[0] += cells
        simplex_totals[1] += corners
    vdw_totals = [0, 0]
    for path in sorted((cert / "vdw").glob("*.json")):
        vdw_totals[0] += 1
        vdw_totals[1] += verify_vdw(path)

    expected = manifest["expected"]
    assert algebraic[0] == expected["algebraic_records"]
    assert algebraic[1] == expected["article_table_records"]
    assert algebraic[2] == expected["algebraic_translates"]
    assert len(list((cert / "simplex").glob("*.json"))) == expected["simplex_files"]
    assert simplex_totals[1] == expected["simplex_corners"]
    assert vdw_totals[0] == expected["vdw_certificates"]
    assert vdw_totals[1] == expected["vdw_progressions"]
    print(
        "ALL CERTIFICATES VERIFIED: "
        f"{algebraic[0]} algebraic, {len(list((cert / 'simplex').glob('*.json')))} simplex, "
        f"{vdw_totals[0]} van der Waerden; zero violations."
    )


if __name__ == "__main__":
    main()
