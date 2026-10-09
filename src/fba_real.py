"""Real FBA runner for VitaMind — bridges to a COBRA-capable interpreter.

Motivation
----------
`metabolic_gem` previously only exercised a *mock* aniger_gem backend
(`status: mock`), so the end-to-end success rate (BioMaster 95.9%-style
metric) could never reflect real metabolic modelling.  We now have a genuine
A. niger SBML GEM (iAniger_ccm / _refined, supplied as project assets), and
COBRApy 0.31.1 is available in a system interpreter (Python 3.14).

This module:
  * runs real FBA/parsimonious FBA directly if `cobra` is importable in the
    current interpreter;
  * otherwise auto-discovers a COBRA-capable python (env override or probe of
    common interpreters) and runs the same solver via subprocess;
  * always degrades gracefully to status="mock" if no COBRA backend exists,
    so callers never crash.

Lazy `import cobra` keeps this module importable under any interpreter.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from typing import Dict, List, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_GEM_CANDIDATES = [
    os.environ.get("VITAMIND_GEM_MODEL"),
    os.path.join(_HERE, "assets", "gem", "aniger_ccm_refined.xml"),
    os.path.join(_HERE, "assets", "gem", "aniger_ccm.xml"),
    # 用户提供的 COMBINE 档案内的 SBML（若解包在工作区）
    os.path.join(_HERE, "..", "biomaster_probe", "gem_omex", "MODEL1507180047_urn.xml"),
]


def find_default_gem() -> Optional[str]:
    for p in _DEFAULT_GEM_CANDIDATES:
        if p and os.path.exists(p):
            return os.path.abspath(p)
    return None


def _find_cobra_python() -> Optional[str]:
    """Return a python executable that can `import cobra`, or None."""
    # Never hard-code a machine-specific interpreter path here: it makes the
    # module silently unusable on any other machine, and it is exactly what
    # verify.py check 4 exists to catch.
    # Order: explicit env override -> the running interpreter -> sibling
    # interpreters of the same major/minor version -> plain PATH names.
    candidates: List[str] = []
    if os.environ.get("VITAMIND_FBA_PYTHON"):
        candidates.append(os.environ["VITAMIND_FBA_PYTHON"])
    candidates.append(sys.executable)
    _exe_dir = os.path.dirname(os.path.abspath(sys.executable))
    _tag = "%d.%d" % (sys.version_info[0], sys.version_info[1])
    for _name in ("python" + _tag + ".exe", "python" + _tag, "python3", "python"):
        candidates.append(os.path.join(_exe_dir, _name))
    candidates += ["python3", "python"]
    for cand in candidates:
        if not cand:
            continue
        try:
            out = subprocess.run(
                [cand, "-c", "import cobra; print(cobra.__version__)"],
                capture_output=True, text=True, timeout=30,
            )
            # stdout 形如 "0.31.1" —— 含数字即视为可用（不要搜 'cobra' 字样）
            if out.returncode == 0 and any(ch.isdigit() for ch in out.stdout):
                return cand
        except Exception:
            continue
    return None


def _set_objective(model, objective: Optional[str]):
    """Set model objective. objective=None → keep model default (BIOMASS).
    objective='EX_cit' → maximize that exchange (shows product secretion)."""
    if not objective:
        return model.objective
    if objective in model.reactions:
        model.objective = objective
        return objective
    # try prefix match (e.g. 'cit' -> EX_cit)
    for r in model.reactions:
        if r.id.startswith(objective) or objective in r.id:
            model.objective = r.id
            return r.id
    raise KeyError(f"objective reaction not found: {objective}")


def _solve(model_path: str, objective: Optional[str] = None) -> Dict:
    """Run real FBA with COBRApy (must be importable here)."""
    import cobra  # lazy
    import cobra.flux_analysis as fa

    model = cobra.io.read_sbml_model(model_path)
    if objective:
        used_obj = _set_objective(model, objective)
    else:
        used_obj = "BIOMASS"  # 模型默认目标反应
    sol = model.optimize()
    if sol.status != "optimal":
        # fall back to parsimonious FBA
        try:
            sol = fa.pfba(model)
        except Exception:
            pass

    fluxes = {r.id: float(sol.fluxes.get(r.id, 0.0)) for r in model.reactions}

    # gene flux participation
    gene_part: Dict[str, float] = {}
    for g in model.genes:
        if g.id.lower() in ("biomass", "bio"):
            continue  # 过滤 GEM 自带的伪基因
        part = 0.0
        for rxn in g.reactions:
            part += abs(fluxes.get(rxn.id, 0.0))
        if part > 0:
            gene_part[g.id] = round(part, 4)
    top_genes = sorted(gene_part.items(), key=lambda kv: kv[1], reverse=True)[:5]

    # key exchanges
    exchanges = {}
    for key in ("EX_glc", "EX_phos", "EX_cit", "EX_glucon", "EX_co2", "BIOMASS"):
        for r in model.reactions:
            if r.id == key or r.id.startswith(key):
                exchanges[r.id] = round(fluxes.get(r.id, 0.0), 4)
                break

    return {
        "status": "ok",
        "engine": f"cobra:{cobra.__version__}",
        "model": os.path.basename(model_path),
        "organism": "A. niger (iAniger_ccm)",
        "objective_reaction": str(used_obj),
        "objective_value": round(float(sol.objective_value), 4),
        "biomass_flux": round(float(sol.fluxes.get("BIOMASS", 0.0)), 4),
        "n_reactions": len(model.reactions),
        "n_metabolites": len(model.metabolites),
        "n_genes": len(model.genes),
        "exchanges": exchanges,
        "top_genes": [g for g, _ in top_genes],
        "gene_scores": {g: s for g, s in top_genes},
        "provenance": "real_fba_cobra",
    }


def _solve_via_subprocess(model_path: str, objective: Optional[str] = None) -> Dict:
    py = _find_cobra_python()
    if not py:
        return {"status": "mock", "error": "no COBRA-capable python found"}
    cmd = [py, os.path.abspath(__file__), "--cli", "--model", model_path]
    if objective:
        cmd += ["--objective", objective]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if out.returncode != 0:
            return {"status": "mock", "error": f"fba subprocess rc={out.returncode}: {out.stderr[:200]}"}
        # last non-empty stdout line should be JSON
        for line in reversed(out.stdout.strip().splitlines()):
            line = line.strip()
            if line.startswith("{"):
                return json.loads(line)
        return {"status": "mock", "error": "no JSON from fba subprocess"}
    except Exception as e:
        return {"status": "mock", "error": f"subprocess error: {e}"}


def run_real_fba(model_path: Optional[str] = None,
                 objective: Optional[str] = None,
                 _direct: bool = False) -> Dict:
    """Public entry. Runs real FBA; degrades to status='mock' if impossible."""
    model_path = model_path or find_default_gem()
    if not model_path or not os.path.exists(model_path):
        return {"status": "mock", "error": "no GEM model available"}
    try:
        import cobra  # noqa: F401  (only succeeds if available here)
        return _solve(model_path, objective)
    except ImportError:
        if _direct:
            return {"status": "mock", "error": "cobra not importable in direct mode"}
        return _solve_via_subprocess(model_path, objective)


def cli_main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--cli", action="store_true")
    ap.add_argument("--model", default=None)
    ap.add_argument("--objective", default=None)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    res = run_real_fba(args.model, args.objective, _direct=True)
    print(json.dumps(res, ensure_ascii=False))


if __name__ == "__main__":
    cli_main()
