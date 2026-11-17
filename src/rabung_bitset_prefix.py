#!/usr/bin/env python3
"""Greedy/random bitset-prefix search for a generalised Rabung set."""
from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path


def factors(n):
    out, q = [], 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out


def primitive_root(p):
    fs = factors(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fs):
            return g
    raise ValueError("no primitive root; p must be prime")


def bits(x):
    while x:
        b = x & -x
        yield b.bit_length() - 1
        x -= b


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("p", type=int)
    ap.add_argument("r", type=int)
    ap.add_argument("t", type=int)
    ap.add_argument("seconds", type=float)
    ap.add_argument("seed", type=int)
    ap.add_argument("output")
    ap.add_argument("--candidates", type=int, default=16)
    args = ap.parse_args()
    p, r, t = args.p, args.r, args.t
    if t < 3 or r < 2 or (p - 1) % r:
        ap.error("require t>=3, r>=2, and r | p-1")
    rng = random.Random(args.seed)
    g = primitive_root(p)
    chi = [r] * p
    x = 1
    for i in range(p - 1):
        chi[x] = i % r
        x = x * g % p
    mask = (1 << p) - 1
    B = [0] * r
    for x in range(1, p):
        B[chi[x]] |= 1 << x

    def rotate(b, s):
        s %= p
        return ((b >> s) | (b << (p - s))) & mask if s else b

    def direct(S):
        bad = 0
        for y in range(p):
            values = [chi[(y + s) % p] for s in S if (y + s) % p]
            if len(set(values)) == 1:
                bad += 1
        return bad

    def good_last(S, I):
        forbidden = 0
        for c in range(r):
            for y in bits(I[c]):
                forbidden |= rotate(B[c], y)  # (A): no zero level
                forbidden |= 1 << ((-y) % p)  # (B): new level is zero
            for s0 in S:
                if all(chi[(s - s0) % p] == c for s in S if s != s0):
                    forbidden |= rotate(B[c], (-s0) % p)  # (C)
        for s in S:
            forbidden |= 1 << s
        return (~forbidden) & mask

    start = time.perf_counter()
    tries = 0
    best_score = 10**9
    last = start
    while time.perf_counter() - start < args.seconds:
        tries += 1
        S = [0, 1]
        inS = {0, 1}
        I = [rotate(B[c], 0) & rotate(B[c], 1) for c in range(r)]
        while len(S) < t - 2:
            chosen = chosenI = score_best = None
            for _ in range(args.candidates):
                z = rng.randrange(2, p)
                if z in inS:
                    continue
                J = [I[c] & rotate(B[c], z) for c in range(r)]
                score = sum(q.bit_count() for q in J)
                if score_best is None or score < score_best:
                    chosen, chosenI, score_best = z, J, score
            if chosen is None:
                break
            S.append(chosen)
            inS.add(chosen)
            I = chosenI
        if len(S) != t - 2:
            continue
        best_score = min(best_score, sum(q.bit_count() for q in I))
        # Exhaust the penultimate element, then use the bitset last-element test.
        for z in range(2, p):
            if z in inS:
                continue
            Sz = S + [z]
            Iz = [I[c] & rotate(B[c], z) for c in range(r)]
            good = good_last(Sz, Iz)
            if good:
                T = Sz + [next(bits(good))]
                n_bad = direct(T)
                if n_bad:
                    raise RuntimeError(f"bitset/direct mismatch: {T}, bad={n_bad}")
                data = {
                    "p": p, "r": r, "t": t, "g": g, "bad": 0, "weights": T,
                    "tries": tries, "seconds": time.perf_counter() - start,
                    "engine": "greedy-bitset-prefix-exhaustive-last-two",
                }
                Path(args.output).write_text(json.dumps(data, indent=2) + "\n")
                print("FOUND", json.dumps(data), flush=True)
                return 0
        if time.perf_counter() - last > 5:
            print(f"PROGRESS p={p} tries={tries} min_prefix_popcount={best_score}", flush=True)
            last = time.perf_counter()
    print(f"DONE p={p} tries={tries} min_prefix_popcount={best_score}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
