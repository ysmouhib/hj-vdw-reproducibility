[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24297573.svg)](https://doi.org/10.5281/zenodo.24297573)
[![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b.svg)](https://arxiv.org/pdf/2611.22205)



# Hales–Jewett and van der Waerden certificates

Reproducible certificates and search code for **More lower bounds for Hales–Jewett and van der Waerden numbers** (Younes Mouhib, 2026).

The strongest retained results include:

- `HJ(5,2) >= 84`, `HJ(3,4) >= 51`, `HJ(4,3) >= 164`, `HJ(5,4) >= 7418`, and `HJ(5,6) >= 34034`;
- `kappa_ow(5,4) >= 7416` and `kappa_ow(5,6) >= 34032`;
- `W(8,3) > 522`, `W(11,3) > 1193`, and `W(17,3) > 3620`.

The claims about `kappa_ow` concern the **one-weight bracket ceiling**, not unrestricted colourings or the Collapse Conjecture.

## Verify everything

Python 3 is sufficient; no solver or third-party package is needed:

```sh
python3 src/verify_all.py
```

A successful run ends with:

```text
ALL CERTIFICATES VERIFIED: 22 algebraic, 3 simplex, 5 van der Waerden; zero violations.
```

The verifier checks manifest hashes and then independently enumerates every algebraic translate, simplex corner, and three-term arithmetic progression. It also reconstructs and hashes every period-`2p` algebraic palette.

## Contents

- `certificates/algebraic.json` — 22 compact `(t,r,p,g,S)` records: 19 article rows, two valid rows omitted in favour of stronger external bounds, and one explicitly superseded `(5,6)` row.
- `certificates/simplex/` — the `T_50` witness, its `T_49` warm start, and the `T_21` calibration witness.
- `certificates/vdw/` — the three article certificates and two clearly marked superseded colourings.
- `certificates/manifest.json` — SHA-256 integrity hashes and expected verification totals.
- `src/verify_all.py` — standalone standard-library verifier.
- `src/` — algebraic, simplex, interval, and exact SAT search code.
- `methods/METHODS.md` — exact algorithms, parameters, commands, and scope of SAT results.
- `REFERENCES.md` — articles, software repositories, and prior-work links.
- `LICENSE` — MIT License, matching the previous project.

No article source, solver logs, or unrelated research files are included.

## Search code

Compile the four C++ programs:

```sh
sh src/build.sh
```

This creates local binaries under `build/`; they are not part of the source release. Exact SAT experiments require:

```sh
python3 -m pip install -r requirements.txt
```

`pycryptosat==5.18.0` supplies the CryptoMiniSat Python interface. `python-sat` is needed only for the original-style global interval CDCL driver. Neither is needed to verify the certificates.

See [methods/METHODS.md](methods/METHODS.md) for reproduction commands, including the deterministic scan that finds the `(5,6)` set `{0,1,654,7769,23051}` at `p=34033`.

## Interpreting search outcomes

A valid JSON certificate is checked independently before use. A time- or conflict-limited SAT result is **UNKNOWN**, never UNSAT. UNSAT inside a Hamming ball excludes only that neighbourhood; simplex UNSAT would concern only symmetric colourings. No next-length timeout is an upper-bound proof.

## Citation and archiving

Cite the article and this repository using `CITATION.cff`. The previous related project is archived at [Zenodo DOI 10.5281/zenodo.21537624](https://doi.org/10.5281/zenodo.21537624). For this repository, create a GitHub release and archive it through Zenodo to obtain a new DOI; do not reuse the previous project's DOI.

## Licence

MIT License, Copyright (c) 2026 Younes Mouhib. This covers the code and certificate collection. References and linked external works retain their own licences.
