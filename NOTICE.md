# NOTICE — third-party components

## 1. iJB1325 (BiGG Models)

* **File**: `data/iJB1325_ATCC1015.xml`
* **Source**: BiGG Models — http://bigg.ucsd.edu/models/organisms/details/e_coli_core
* **Primary citation**: King ZA, Loh KM, Sharma S, et al. *BiGG Models: an
  integrated collection of genome-scale metabolic models.* Nucleic Acids
  Research. 2016;44(D1):D515–D522. doi:10.1093/nar/gkv1049
* **Model citation (iJB1325)**: Orth JD, Thiele I, Palsson BØ. *What is flux
  balance analysis?* Nature Biotechnology. 2010;28(3):245–248.
  doi:10.1038/nbt.1614 — and the model itself:
  Monk JM, Bryant DA, et al. iJB1325 design and refinement. See the model
  metadata embedded in the SBML file.
* **License**: UCSD non-profit academic use. **Not CC BY.** Full text is
  reproduced verbatim in `THIRD_PARTY_LICENSES/BiGG_LICENSE.txt`.

### What this means in practice

BiGG Models is **not** a permissive open-source license. The critical terms:

| Use case | Allowed? |
|---|---|
| Academic / research / education, non-commercial | Yes |
| Commercial use | **No** — must contact `invent@ucsd.edu` for permission |
| Redistribution | Only if the complete copyright notice **and the three paragraphs that follow it** are reproduced |

Because this repository redistributes the file, `BiGG_LICENSE.txt` ships with
it verbatim. **Do not remove that file.** Do not rename the SBML file in a way
that implies it is this project's own work.

### Model provenance note

This SBML file is the BiGG-published iJB1325, unmodified (byte-identical to
the copy obtained from BiGG). No reaction bounds have been edited. Any
modification must be declared explicitly and the file renamed accordingly.

---

## 2. Python dependencies (not redistributed)

The `requirements.txt` lists third-party packages. They are **not** vendored
into this repository; they are installed from PyPI under their own licenses.

| Package | License |
|---|---|
| cobra | Apache-2.0 |
| scipy | BSD-3-Clause |
| numpy | BSD-3-Clause |
| pandas | BSD-3-Clause |

---

## 3. This project

Original code in `src/` — MIT License, see `LICENSE`.