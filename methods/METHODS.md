# Methods and reproduction notes

This file documents the algorithms precisely. The only requirement for checking the mathematical certificates is Python 3; `src/verify_all.py` uses no third-party package.

## 1. Generalised Rabung certificates

Let `p` be prime, `r | (p-1)`, and `g` a primitive root. For nonzero `x`, define `chi(x) = ind_g(x) mod r`. A `t`-element set `S` satisfies condition (R) if, for every `y` in `F_p`, the nonzero values `y+s`, `s in S`, have at least two `chi`-colours.

The article uses the period-`p` colouring for lines with fewer than `p` active coordinates. With integer weights congruent to `S` modulo `p`, it gives `HJ(t,r) >= p` and `kappa_ow(t,r) >= p-1`. If the weights have both parities, the period-`2p` palette

- `chi(x mod p)` when `p` does not divide `x`;
- colour `0` at `x=0` and colour `1` at `x=p`;

proves `HJ(t,r) >= p+1` in dimension `p`. Every algebraic JSON record contains the complete compact certificate `(t,r,p,g,S)` and the SHA-256 hash of this full palette.

### Bitset-prefix elimination

For each colour `c`, let `B_c` be the bitset of residues with `chi(x)=c`. For a prefix `S'`, form

`I_c = intersection_{s in S'} (B_c-s)`.

A final element `x` is illegal exactly when at least one of these mechanisms occurs:

- **(A)** `x in B_c-y` for some `y in I_c`: no level vanishes and the translate is monochromatic.
- **(B)** `x=-y` for some `y in I_c`: the new level is the exceptional zero.
- **(C)** `x in B_c+s0`, where `s0` is a prefix element and all other prefix differences have colour `c`: the level of `s0` is zero.

Affine normalisation permits prefixes containing `{0,1}`. `rabsearch_exhaustive.cpp` is the actual five-element deterministic scanner: it enumerates `{0,1,a,b}` and finds the least legal `x>b`. A hit is checked again by direct enumeration of all `p` translates before it is written.

Reproduce the strongest `(5,6)` certificate (the scan stops at its first hit):

```sh
g++ -O3 -std=c++17 src/rabsearch_exhaustive.cpp -o /tmp/rabsearch_exhaustive
/tmp/rabsearch_exhaustive 34033 6 /tmp/reproduced_5_6.json
```

Expected standard output:

```text
FOUND p=34033 r=6 S=0,1,654,7769,23051,
```

The consolidation run examined 21,981,797 prefixes and performed 1,214,640,303 bitset union/complement operations before this hit. Runtime is hardware-dependent; the mathematical check is the independent verifier, not the runtime.

`rabung_bitset_prefix.py` implements greedy prefix growth. At each level it samples candidate elements and minimises `sum_c |I_c|`, then exhausts the penultimate element and applies (A)–(C) to the last one. Known-answer control:

```sh
python3 src/rabung_bitset_prefix.py 163 3 4 10 1 /tmp/p163.json
# Returns S = {0,1,20,21}.
```

`rabung_weight_local.cpp` is the historical weighted local-search implementation on `S` itself. It maintains translate/colour counts incrementally, samples a high-weight bad translate from eight draws, tests 32 replacements of a non-normalised element, and uses breakout weights, perturbations, and restarts. Its acceptance and noise rules are retained verbatim in the source for provenance; only the command-line argument-count check was corrected. Example control starting one element away from a solution:

```sh
g++ -O3 -std=c++17 src/rabung_weight_local.cpp -o /tmp/rabung_weight_local
/tmp/rabung_weight_local 163 3 4 2 1 /tmp/local163.json 0,1,20,22
```

## 2. Symmetric simplex search

For three letters, row `a` and position `b` store the colour of `(a,b,n-a-b)`. The corner triple with scale `k` and base `(a,b)` is

`{(a+k,b), (a,b+k), (a,b)}`.

`simplex_warm_weighted.cpp` maintains all corner/colour counts and a constant-time list of monochromatic corners. Its actual parameters are:

- tabu tenure uniformly from 12 through 40, with aspiration when a move beats the global best;
- eight random bad-corner samples, retaining the one of largest breakout weight;
- a 25% breakout-weight increment when the best weighted gain is nonnegative;
- a 1% random move when the best gain is positive;
- random tie-breaking, perturbation on restart, and a restart after 2,000,000 steps without global improvement.

There is **no periodic weight-smoothing operation** in this implementation.

A warm start in direction `d` assigns a cell `v` of `T_n` the old colour of `v-e_d` at the first available cyclic direction `d,d+1,d+2`. This preserves all corners lying inside the copied face. Run from the supplied `T_49` certificate with an explicit new seed:

```sh
g++ -O3 -std=c++17 src/simplex_warm_weighted.cpp -o /tmp/simplex_warm_weighted
/tmp/simplex_warm_weighted 50 60 1 /tmp/T50_new.json \
  certificates/simplex/HJ34_T49_warm_start.json 1
```

This is a stochastic search command, not a promised runtime. The retained `T_50` witness predates this packaging and its original random seed was not logged; its 1,326 digits and all 22,100 corners are nevertheless completely checked by `verify_all.py`. No seed or runtime has been invented for the old run.

## 3. Van der Waerden local search

The original searches began with Heule's unrolled cyclic colourings of lengths 510, 1182, and 3546 for 8, 11, and 17 colours. Endpoint extension therefore first produced lengths 511, 1183, and 3547. At each extension, one endpoint integer receives a random colour and the **whole** interval is repaired; no old position is fixed.

`vdw_tabu.cpp` is a clean implementation of this method for three-term APs. It:

1. selects a monochromatic AP;
2. evaluates recolouring each of its three positions with every other colour;
3. minimises the incremental change in the number of monochromatic APs, with random tie-breaking;
4. keeps a changed position tabu for 2 through 8 steps, with aspiration;
5. makes 2% random moves;
6. perturbs the best state after 200,000 steps without a global improvement.

Example with a supplied warm colouring:

```sh
g++ -O3 -std=c++17 src/vdw_tabu.cpp -o /tmp/vdw_tabu
/tmp/vdw_tabu 523 8 60 1 /tmp/W83_N523.json \
  certificates/vdw/W_8_3_N522.json
```

Add a final `1` to prepend instead of append when the warm file has length `N-1`. A time-capped exit without `bad=0` is UNKNOWN. The final certificates are valid independently of whether a fresh stochastic run finds them again.

## 4. Exact SAT methods

### 4.1 CryptoMiniSat simplex encoding

The original exact environment used CryptoMiniSat/`pycryptosat` 5.18.0, one solver thread, Python 3, and C++ programs compiled with optimisation. For cell `i` and colour `c`, introduce `x_(i,c)`. Add:

- one at-least-one-colour clause per cell;
- for every corner tuple and colour, one clause forbidding all its cells from having that colour.

At-most-one clauses are unnecessary for witness extraction: choosing one true colour per cell cannot create a monochromatic corner. The original `T_21`, three-colour calibration used one additional unit clause fixing cell 0 to colour 0:

```sh
python3 -m pip install pycryptosat==5.18.0
python3 src/simplex_cms.py --t 3 --r 3 --n 21 --fix-first --time-limit 600
```

The optional `--symmetry` flag adds complete first-occurrence breaking:

`x_(i,c) -> OR_(h<i) x_(h,c-1)` for `c>=1`.

This stronger scheme is available for new experiments; it was not part of the original unit-clause calibration.

### 4.2 One-weight encoding

For weights `omega`, only the realised levels `ell = <omega,type>` need variables. In the two-colour case, one Boolean per level and two clauses per distinct line homothet suffice. The released encoder reproduces the 64-level, 379-homothet calibration instance:

```sh
python3 src/oneweight_cms.py --omega 0,2,3,5 --r 2 --n 13 \
  --time-limit 120 --output /tmp/HJ42_oneweight.json
```

For `S={0,1,3,11}`, the nonperiodic weights `(0,2,3,5)` are `S+2` modulo 13; the period-26 representative instead uses `(0,1,3,11)` with the exceptional palette values specified in Section 1. Both representations certify `HJ(4,2)>=14`.

### 4.3 Exact Hamming-neighbourhood repair

`simplex_cms.py --warm W --radius D` adds at-most-one-colour clauses and a sequential-counter cardinality constraint on the literals `not x_(i,warm_i)`. Thus the selected ordinary colouring differs from the warm one in at most `D` cells. Do not combine a warm-centred Hamming ball with colour-renaming symmetry.

A completed UNSAT call here proves only that the specified Hamming ball contains no solution. It does not prove simplex UNSAT and says nothing about unrestricted Hales--Jewett colourings. For example, this command checks the zero-radius control around the supplied certificate:

```sh
python3 src/simplex_cms.py --t 3 --r 4 --n 50 \
  --warm certificates/simplex/HJ34_T50.json --radius 0 --time-limit 60
```

### 4.4 Global interval SAT and budgets

`vdw_cms.py` supplies the analogous unrestricted one-hot encoding for CryptoMiniSat. `vdw_global_pysat.py` retains the original-style CDCL driver with warm colourings as **phase hints only**, no fixed positions, and optional first-occurrence symmetry. The strongest retained global runs did not use colour symmetry.

Conflict and time limits are explicit. `SAT` is followed by direct certificate checking; completed `UNSAT` refers only to the encoded instance; `None` from a budgeted call is always reported as `UNKNOWN`, never as `UNSAT`.

## 5. Verification summary

Run:

```sh
python3 src/verify_all.py
```

The verifier checks file integrity, primality, divisibility, primitive roots, set cardinality and mixed parity, every translate in condition (R), every period-`2p` palette hash, every simplex corner, and every three-term AP. The full collection comprises 22 algebraic records (19 article rows, two valid rows omitted because external bounds are stronger, and one superseded `(5,6)` row), three simplex colourings, and five interval colourings. Only records marked `article_table` or `article_result` support the main result tables; superseded records are retained solely for audit history.
