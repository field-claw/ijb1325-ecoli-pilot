# iJB1325 pilot — minimal reproducible Escherichia coli core

A deliberately small pilot repository: one genome-scale metabolic model, the
four modules that actually use it, and a self-verification script that
reproduces the reference numbers or fails loudly.

Scope was chosen to answer one question honestly: **is a public release
feasible at all?** The answer is yes, at 6.4 MB, with one licence obligation
attached.

## Why this exists

The full project is ~1.7 GB across thousands of files and touches six genome
models with unclear redistribution terms. That is a bad first upload. This
repository is the minimal slice that a stranger can clone and check in under
ten seconds.

## What is here

```
data/iJB1325_ATCC1015.xml   BiGG iJB1325, unmodified (6.4 MB)
src/fba_real.py             real LP FBA via COBRApy
src/enzyme_constrained_gem.py   enzyme-constrained layer
src/ecoli_darkmatter_gap.py     kcat/EC annotation coverage analysis
src/ecoli_darkmatter_eval.py    evaluation for the above
verify.py                   self-verification, fail-closed
```

## Quick start

```bash
python -m pip install -r requirements.txt
python verify.py
```

Expected output (measured on this machine, cobra 0.31.1, GLPK):

```
  PASS license file present               BiGG_LICENSE.txt 3999 bytes
  PASS src/ free of absolute paths        4 files scanned, 0 hits
  PASS model dimensions                   2320rxn/1818met/1325gene
  PASS pFBA growth rate                   0.9399 vs 0.9399 1/h
  PASS carbon guardrail                   growth 0.9399 -> 0.0000 (100.0% drop)

RESULT: PASS -- all 5 checks passed.
```

If cobra or scipy is missing, checks 1–3 report `BLOCKED` rather than `PASS`.
An unrunnable check is never reported as a success.

## What the numbers mean

| Number | Value | What it is |
|---|---|---|
| Reactions / metabolites / genes | 2320 / 1818 / 1325 | Published iJB1325 dimensions — this file *is* iJB1325 |
| Growth rate | 0.9399 h⁻¹ | pFBA optimum on glucose minimal medium |
| Carbon guardrail | 0.9399 → 0.0000 (100% drop) | Closing all 353 carbon-bearing boundaries kills growth |

The guardrail is the one that matters. A model that keeps growing after its
carbon supply is cut off has an internal carbon loop, and every flux
prediction from it is fiction. This one does not.

## What the numbers do NOT mean

Stated plainly, because the failure mode here is not being wrong, it is
being *narrowly* right and letting someone assume more:

- **0.9399 h⁻¹ is one machine, one solver, one medium.** It is not an
  experimental measurement and it is not a claim about real E. coli.
- **Matching published dimensions proves file identity, not scientific
  validity.** A model can be byte-identical to BiGG and still be wrong.
- **The carbon guardrail proves the flux accounting is carbon-bounded.** It
  predicts no phenotype whatsoever.
- **No gene essentiality, no dFBA, no multi-strain community results are
  included here**, despite existing internally. Those are out of scope for
  the pilot on purpose.

## Licence — read this before commercial use

- `src/` — MIT, see `LICENSE`.
- `data/iJB1325_ATCC1015.xml` — **BiGG Models licence, which is NOT CC BY.**
  Academic and non-commercial use only. Commercial use requires prior
  permission from `invent@ucsd.edu`. Full text in
  `THIRD_PARTY_LICENSES/BiGG_LICENSE.txt`; see `NOTICE.md`.

The BiGG licence permits redistribution only if the complete copyright notice
**and the three paragraphs that follow it** are reproduced. That file ships
here verbatim. Do not remove it, and do not present the model as this
project's own work.

## A portability trap worth knowing

iJB1325 prefixes its exchange reactions with `BOUNDARY_`, not `EX_`. Code
that scans for `EX_*` gets an empty list, silently skips its carbon guardrail,
and reports a pass. That exact bug was measured while building this
repository — the guardrail reported "0% drop" until the prefix was fixed.
`verify.py` fails loudly when the boundary count looks wrong for this reason.

## Citations

- King ZA, Loh KM, Sharma S, et al. *BiGG Models: an integrated collection of
  genome-scale metabolic models.* Nucleic Acids Res. 2016;44(D1):D515–D522.
  doi:10.1093/nar/gkv1049
- Orth JD, Thiele I, Palsson BØ. *What is flux balance analysis?* Nat
  Biotechnol. 2010;28(3):245–248. doi:10.1038/nbt.1614
- Monod J. *The growth of bacterial cultures and the relation between the
  specific growth rate and the concentration of limiting substrate.*
  Ann Inst Pasteur. 1949;37:8–111.

## Scope of the pilot

This repository answers "can we ship it". It does not attempt to ship
everything. Models whose redistribution terms could not be verified were
deliberately excluded rather than included with a hedge.