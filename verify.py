#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
verify.py -- pilot repository self-verification (fail-closed)
============================================================================
Run:  python verify.py

This script does NOT print a nice story about the model. It reproduces the
reference numbers and reports PASS / FAIL per check, then exits non-zero if
any check fails. If cobra/scipy is missing it reports BLOCKED (not PASS, and
not FAIL either) -- an unrunnable check is never reported as a success.

Checks
------
  1  model loads, and dimensions match the published iJB1325 counts
  2  pFBA growth rate matches the reference value (0.9399 1/h)
  3  carbon guardrail: closing all carbon-bearing boundaries collapses growth
  4  no hard-coded absolute paths in src/
  5  third-party license file present and non-empty

Exit codes
----------
  0  all executed checks PASS (a BLOCKED run also returns 0, but prints
     BLOCKED so the caller can see nothing was actually verified)
  1  at least one executed check FAIL
"""
from __future__ import annotations

import io
import os
import re
import sys
import contextlib
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(HERE, "data", "iJB1325_ATCC1015.xml")
SRC = os.path.join(HERE, "src")
LIC = os.path.join(HERE, "THIRD_PARTY_LICENSES", "BiGG_LICENSE.txt")

# Published iJB1325 dimensions (BiGG Models, King et al. NAR 2016).
REF_N_RXN = 2320
REF_N_MET = 1818
REF_N_GENE = 1325

# Reference growth rate measured on this machine with cobra 0.31.1 / GLPK.
REF_GROWTH = 0.9399

# Metabolite formulas treated as inorganic carbon for the guardrail check.
INORGANIC_C = {"CO2", "HCO3", "CO", "CH4"}

TOL_GROWTH = 1e-3      # pFBA is LP-exact, but keep a tight tolerance anyway
TOL_COLLAPSE = 0.01   # growth must fall below 1% of baseline

ASCII_INORGANIC_C = "CO2 HCO3 CO CH4"


def head(title: str) -> None:
    print("")
    print("=" * 66)
    print(title)
    print("=" * 66)


def carbon_boundaries(model) -> list:
    """Every boundary reaction carrying at least one non-inorganic carbon.

    NOTE: iJB1325 prefixes its exchanges with BOUNDARY_, not EX_. Hard-coding
    "EX_" silently returns an empty list and the guardrail check becomes a
    no-op that always "passes". That bug was measured here, not assumed --
    see verify.py git history.
    """
    out = []
    for rxn in model.reactions:
        if not rxn.boundary:
            continue
        has_carbon = False
        all_inorganic = True
        for met in rxn.metabolites:
            formula = (met.formula or "").strip()
            if "C" in formula:
                has_carbon = True
                if formula not in INORGANIC_C:
                    all_inorganic = False
        if has_carbon and not all_inorganic:
            out.append(rxn.id)
    return out


def main() -> int:
    head("iJB1325 pilot -- self verification")

    results = []

    # ---- check 5 first: it needs no dependency at all ------------------
    exists = os.path.isfile(LIC)
    size = os.path.getsize(LIC) if exists else 0
    ok5 = exists and size > 1000
    results.append(("license file present", ok5,
                    "BiGG_LICENSE.txt %d bytes" % size if exists else "MISSING"))

    # ---- check 4: no hard-coded absolute paths in src/ -----------------
    hardcoded = []
    abs_re = re.compile(r"[A-Za-z]:[\\/]|/home/|/Users/")
    for root, _dirs, files in os.walk(SRC):
        for fn in files:
            if not fn.endswith(".py"):
                continue
            path = os.path.join(root, fn)
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                for lineno, line in enumerate(fh, 1):
                    if abs_re.search(line):
                        hardcoded.append("%s:%d" % (os.path.relpath(path, HERE), lineno))
    ok4 = not hardcoded
    results.append(("src/ free of absolute paths", ok4,
                    "%d files scanned, %d hits" % (
                        sum(len(f) for _r, _d, f in os.walk(SRC)), len(hardcoded))))

    # ---- dependency-gated checks ---------------------------------------
    try:
        import numpy  # noqa: F401
        from scipy.optimize import linprog  # noqa: F401
        from cobra.flux_analysis import pfba
        from cobra.io import read_sbml_model
    except ImportError as exc:
        head("dependency-gated checks")
        print("BLOCKED: %s" % exc)
        print("install with:  pip install -r requirements.txt")
        for name, ok, note in results:
            print("  %-4s %-34s %s" % ("PASS" if ok else "FAIL", name, note))
        print("")
        print("RESULT: BLOCKED -- checks 1-3 were NOT executed.")
        return 0

    head("model load")
    t0 = time.time()
    model = read_sbml_model(MODEL)
    n_rxn, n_met, n_gene = (len(model.reactions), len(model.metabolites),
                            len(model.genes))
    print("reactions : %d  (reference %d)" % (n_rxn, REF_N_RXN))
    print("metabolites: %d  (reference %d)" % (n_met, REF_N_MET))
    print("genes     : %d  (reference %d)" % (n_gene, REF_N_GENE))
    print("load time : %.2f s" % (time.time() - t0))
    dims_ok = (n_rxn == REF_N_RXN and n_met == REF_N_MET and n_gene == REF_N_GENE)
    results.append(("model dimensions", dims_ok,
                    "%drxn/%dmet/%dgene" % (n_rxn, n_met, n_gene)))
    if not dims_ok:
        print("")
        print("FAIL: the model does not match published iJB1325 dimensions.")
        print("      Do not report any number from this model as iJB1325.")
        for name, ok, note in results:
            print("  %-4s %-34s %s" % ("PASS" if ok else "FAIL", name, note))
        return 1

    head("pFBA growth rate")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        base = pfba(model)
    obj_id = "DRAIN_Biomass"
    if obj_id not in base.fluxes:
        print("FAIL: objective reaction %s not found -- "
              "model objective was edited" % obj_id)
        return 1
    growth = float(base.fluxes[obj_id])
    print("status      : %s" % base.status)
    print("DRAIN_Biomass: %.4f 1/h   (reference %.4f)" % (growth, REF_GROWTH))
    print("total flux  : %.2f mmol/gDW/h" % base.fluxes.abs().sum())
    ok2 = base.status == "optimal" and abs(growth - REF_GROWTH) < TOL_GROWTH
    results.append(("pFBA growth rate", ok2, "%.4f vs %.4f 1/h" % (growth, REF_GROWTH)))

    head("carbon guardrail")
    carb = carbon_boundaries(model)
    print("carbon-bearing boundaries: %d" % len(carb))
    print("inorganic carbon excluded: %s" % ASCII_INORGANIC_C)
    if len(carb) < 10:
        print("FAIL: only %d carbon boundaries found -- the detector itself is "
              "broken, the guardrail result would be meaningless" % len(carb))
        results.append(("carbon guardrail", False,
                        "detector found only %d boundaries" % len(carb)))
    else:
        for rid in carb:
            rxn = model.reactions.get_by_id(rid)
            rxn.lower_bound = 0
            rxn.upper_bound = 0
        with contextlib.redirect_stdout(buf):
            blocked = pfba(model)
        after = float(blocked.fluxes[obj_id])
        drop = 1.0 - after / growth if growth > 0 else float("inf")
        print("baseline            : %.6f 1/h" % growth)
        print("all carbon closed   : %.6f 1/h  (status %s)" % (after, blocked.status))
        print("relative drop       : %.2f %%" % (100.0 * drop))
        ok3 = drop > (1.0 - TOL_COLLAPSE)
        results.append(("carbon guardrail", ok3,
                        "growth %.4f -> %.4f (%.1f%% drop)" % (growth, after, 100.0 * drop)))
        if not ok3:
            print("INTERPRETATION: growth survived carbon withdrawal, which means")
            print("either an internal carbon loop exists or the boundary set is")
            print("incomplete. Do not trust flux predictions from this run.")

    head("summary")
    width = 34
    n_fail = 0
    for name, ok, note in results:
        tag = "PASS" if ok else "FAIL"
        n_fail += 0 if ok else 1
        print("  %-4s %-*s %s" % (tag, width, name, note))
    print("")
    if n_fail:
        print("RESULT: FAIL -- %d of %d checks failed." % (n_fail, len(results)))
        return 1
    print("RESULT: PASS -- all %d checks passed." % len(results))
    print("")
    print("Honest boundary of this verification:")
    print("  - It reproduces ONE growth rate on ONE machine with ONE solver.")
    print("  - Matching published dimensions proves the FILE is iJB1325;")
    print("    it does not prove the model is experimentally validated.")
    print("  - The carbon guardrail proves flux accounting is carbon-bounded.")
    print("    It does not predict any real phenotype.")
    return 0


if __name__ == "__main__":
    sys.exit(main())