# -*- coding: utf-8 -*-
"""ecoli_darkmatter_gap —— 量化 E. coli ecGEM 的「功能注释 / kcat 双缺口」。

动因（对比学习：钟博子韬 Bozitao Zhong, Flash Talk "Virtual E. coli", Westlake AIVC Week 2026）：
    钟博子韬的 Venus 生态（VenusX 细粒度功能注释 / MAPred EC 自回归预测 /
    VenusRXN 反应条件酶发现 / VenusMine / VenusFactory2）从**蛋白语言模型 + 结构**
    侧补齐酶功能「暗物质」；VitaMind 栈内 eciML1515 走**真实 kcat 驱动的 GECKO 酶约束**。
    两者能否互补，取决于一个可量化的问题：
        E. coli GEM 里究竟有多少反应是「有基因但无 EC 号 / 无 kcat」？
    这正是 Venus 路线能直接注靶、且能提升 GECKO 覆盖率的地方。

本模块纯标准库（只读 JSON + 解析既有 Markdown 清单），零外部依赖、零网络、确定性。
默认 dry_run：只评估、不修改任何模型。

用法：
    python vitamind_core/ecoli_darkmatter_gap.py [--json]
"""
from __future__ import annotations
import os, re, json, argparse
from typing import Dict, Any

_THIS = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_THIS)
MODEL_PATH = os.path.join(_THIS, "data", "gem", "eciML1515.json")
TOP50_MD = os.path.join(_ROOT, "质谱优先测蛋白清单_eciML1515_Top50_2026-09-07.md")
OUT_JSON = os.path.join(_ROOT, "data", "ecoli_darkmatter_gap.json")

# 剥离既有酶层（pro__L 伪代谢物相关），与 gecko_eciML1515.py 保持一致
_ENZ_LAYER_FRAG = ("pro__L", "PROtex", "PROt2rpp")

# ecModels 同工酶/可逆拆分后缀（需连剥两次：_num1_reverse → _num1 → base）
_RE_SPLIT_SUFFIX = re.compile(r"(_reverse|_rev|_num\d+)$", re.I)


def _is_enzyme_layer(rid: str) -> bool:
    return any(f in rid for f in _ENZ_LAYER_FRAG)


def _has_ec(ann: Dict[str, Any]) -> bool:
    """判定反应是否有 EC 号注释（ec-code 或 rhea/kegg 间接可用作 EC 检索锚点）。"""
    if not isinstance(ann, dict):
        return False
    if ann.get("ec-code"):
        return True
    return False


def _has_rhea_kegg(ann: Dict[str, Any]) -> bool:
    """有 Rhea/KEGG 反应 ID → 可作为 VenusRXN 这类「反应条件酶发现」的输入锚点。"""
    if not isinstance(ann, dict):
        return False
    return bool(ann.get("rhea") or ann.get("kegg.reaction") or ann.get("metanetx.reaction"))


def analyze() -> Dict[str, Any]:
    model = json.load(open(MODEL_PATH, encoding="utf-8"))
    kcat = model.get("enzyme_constraint", {}).get("kcat_MW", {}) or {}
    rxns = [r for r in model.get("reactions", []) if not _is_enzyme_layer(r.get("id", ""))]

    n_total = len(rxns)
    n_gpr = n_kcat = n_ec = 0
    # 2x2x2 交叉：GPR × kcat × EC
    cross: Dict[str, int] = {}
    dark_full, kcat_gap_only, ec_gap_only = [], [], []

    for r in rxns:
        rid = r.get("id", "")
        gpr = (r.get("gene_reaction_rule") or "").strip()
        has_g = bool(gpr)
        has_k = rid in kcat
        ann = r.get("annotation") or {}
        has_e = _has_ec(ann)
        n_gpr += has_g
        n_kcat += has_k
        n_ec += has_e
        key = f"gpr{int(has_g)}_kcat{int(has_k)}_ec{int(has_e)}"
        cross[key] = cross.get(key, 0) + 1
        rec = {"id": rid, "name": r.get("name", ""), "gpr": gpr,
               "has_ec": has_e, "has_rhea_kegg": _has_rhea_kegg(ann)}
        if has_g and not has_e and not has_k:
            dark_full.append(rec)          # 完全暗物质：有基因、无 EC、无 kcat → Venus 双可补
        elif has_g and has_e and not has_k:
            kcat_gap_only.append(rec)      # 有 EC 无 kcat → 需 kcat 预测（非 Venus 直接靶）
        elif has_g and not has_e and has_k:
            ec_gap_only.append(rec)        # 有 kcat 无 EC → 注释缺失但参数已有

    # 孤儿反应：代谢主干上没有任何基因关联（可能是转运/扩散/交换，也可能是真孤儿）
    orphans = [{"id": r.get("id", ""), "name": r.get("name", ""),
                "has_ec": _has_ec(r.get("annotation") or {}),
                "has_rhea_kegg": _has_rhea_kegg(r.get("annotation") or {})}
               for r in rxns if not (r.get("gene_reaction_rule") or "").strip()]

    # Top50 瓶颈酶清单：统计其中「无 EC 号」的条目（这些是酶池占用最大却无功能注释者）
    top50_no_ec, top50_total = [], 0
    if os.path.exists(TOP50_MD):
        for line in open(TOP50_MD, encoding="utf-8"):
            m = re.match(r"^\|\s*(\d+)\s*\|\s*`([^`]+)`\s*\|([^|]*)\|([^|]*)\|([^|]*)\|([^|]*)\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|", line)
            if not m:
                continue
            top50_total += 1
            rank, rid, name, genes, uniprot, ec, pct, cum = m.groups()
            if ec.strip() in ("—", "-", "", "—)"):
                top50_no_ec.append({"rank": int(rank), "id": rid, "name": name.strip(),
                                    "genes": genes.strip(), "uniprot": uniprot.strip(),
                                    "pool_pct": float(pct), "cum_pct": float(cum)})

    # —— 去重口径：合并 _numN / _reverse 同工酶拆分变体 ——
    # 诚实要点：eciML1515 是 ecModels 格式，66% 反应带拆分后缀，直接计数会把同一
    # 生化反应重复计 N 次。此处按 base reaction 归并（任一变体有 kcat/EC 即视为已覆盖）。
    from collections import defaultdict
    b_kcat: Dict[str, bool] = defaultdict(bool)
    b_ec: Dict[str, bool] = defaultdict(bool)
    b_gpr: Dict[str, bool] = defaultdict(bool)
    n_suffixed = 0
    for r in rxns:
        rid = r.get("id", "")
        if _RE_SPLIT_SUFFIX.search(rid):
            n_suffixed += 1
        b = _RE_SPLIT_SUFFIX.sub("", _RE_SPLIT_SUFFIX.sub("", rid))
        ann = r.get("annotation") or {}
        b_kcat[b] = b_kcat[b] or (rid in kcat)
        b_ec[b] = b_ec[b] or _has_ec(ann)
        b_gpr[b] = b_gpr[b] or bool((r.get("gene_reaction_rule") or "").strip())
    n_base = len(b_gpr)
    n_b_enz = sum(1 for b in b_gpr if b_gpr[b])
    n_b_kcat = sum(1 for b in b_gpr if b_kcat[b])
    n_b_ec = sum(1 for b in b_gpr if b_ec[b])
    b_dark = sum(1 for b in b_gpr if b_gpr[b] and not b_ec[b] and not b_kcat[b])
    b_kgap = sum(1 for b in b_gpr if b_gpr[b] and b_ec[b] and not b_kcat[b])
    b_egap = sum(1 for b in b_gpr if b_gpr[b] and not b_ec[b] and b_kcat[b])
    dedup = {
        "note": "按 base reaction 归并 _numN/_reverse 同工酶拆分变体；任一变体有 kcat/EC 即视为覆盖",
        "n_base": n_base, "n_enzymatic": n_b_enz, "n_with_kcat": n_b_kcat, "n_with_ec": n_b_ec,
        "kcat_coverage_of_enzymatic": round(n_b_kcat / n_b_enz, 4) if n_b_enz else None,
        "ec_coverage_of_enzymatic": round(n_b_ec / n_b_enz, 4) if n_b_enz else None,
        "kcat_coverage_of_all": round(n_b_kcat / n_base, 4) if n_base else None,
        "ec_coverage_of_all": round(n_b_ec / n_base, 4) if n_base else None,
        "dark_full": b_dark, "kcat_gap_only": b_kgap, "ec_missing_only": b_egap,
        "dark_frac_of_enzymatic": round(b_dark / n_b_enz, 4) if n_b_enz else None,
        "n_suffixed_variants": n_suffixed,
        "suffix_frac_of_backbone": round(n_suffixed / n_total, 4) if n_total else None,
    }

    return {
        "dedup": dedup,
        "model": {"id": model.get("id"), "n_reactions_raw": len(model.get("reactions", [])),
                  "n_reactions_backbone": n_total, "n_genes": len(model.get("genes", [])),
                  "n_kcat_MW": len(kcat), "enzyme_upperbound": model.get("enzyme_constraint", {}).get("upperbound")},
        "counts": {"n_total": n_total, "n_with_gpr": n_gpr, "n_with_kcat": n_kcat, "n_with_ec": n_ec},
        "coverage": {
            "kcat_coverage_of_all": round(n_kcat / n_total, 4) if n_total else None,
            "kcat_coverage_of_enzymatic": round(n_kcat / n_gpr, 4) if n_gpr else None,
            "ec_coverage_of_enzymatic": round(n_ec / n_gpr, 4) if n_gpr else None,
        },
        "cross_gpr_kcat_ec": cross,
        "dark_matter_full": {"n": len(dark_full), "examples": dark_full[:15]},
        "kcat_gap_only": {"n": len(kcat_gap_only), "examples": kcat_gap_only[:10]},
        "ec_gap_only": {"n": len(ec_gap_only), "examples": ec_gap_only[:10]},
        "orphan_reactions": {"n": len(orphans),
                             "n_with_rhea_kegg": sum(1 for o in orphans if o["has_rhea_kegg"])},
        "top50_bottleneck": {"n_listed": top50_total, "n_without_ec": len(top50_no_ec),
                             "pool_pct_without_ec": round(sum(x["pool_pct"] for x in top50_no_ec), 2),
                             "items": top50_no_ec},
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="E. coli ecGEM 功能注释/kcat 双缺口量化")
    ap.add_argument("--json", action="store_true", help="输出 JSON 到 stdout")
    args = ap.parse_args()

    res = analyze()
    c, cov = res["counts"], res["coverage"]
    print("=== E. coli ecGEM 功能注释 / kcat 双缺口量化（eciML1515 真实资产）===")
    print(f"  主干反应 {c['n_total']} ｜ 基因 {res['model']['n_genes']} ｜ kcat_MW {res['model']['n_kcat_MW']} "
          f"｜ 酶容量上界 {res['model']['enzyme_upperbound']}")
    print(f"  有 GPR（酶促）反应   : {c['n_with_gpr']}")
    print(f"  有 kcat_MW 反应      : {c['n_with_kcat']}  (占主干 {cov['kcat_coverage_of_all']:.2%} / 占酶促 {cov['kcat_coverage_of_enzymatic']:.2%})")
    print(f"  有 EC 号注释反应     : {c['n_with_ec']}  (占酶促 {cov['ec_coverage_of_enzymatic']:.2%})")
    print()
    print("  —— Venus 路线靶心 ——")
    print(f"  ① 完全暗物质（有GPR·无EC·无kcat）: {res['dark_matter_full']['n']} 条  → VenusX/MAPred 补 EC + kcat 预测")
    print(f"  ② 有EC·无kcat（参数缺口）        : {res['kcat_gap_only']['n']} 条  → 需 kcat 预测/DLKcat，非 Venus 直接靶")
    print(f"  ③ 有kcat·无EC（仅注释缺失）      : {res['ec_gap_only']['n']} 条  → VenusX 可回填注释，零成本增益")
    print(f"  ④ 孤儿反应（无GPR）              : {res['orphan_reactions']['n']} 条 "
          f"(其中 {res['orphan_reactions']['n_with_rhea_kegg']} 条有 Rhea/KEGG 锚点 → VenusRXN 可招募酶)")
    d = res["dedup"]
    print()
    print("  —— 口径校正：eciML1515 是 ecModels 格式，66% 反应带 _numN/_reverse 同工酶拆分后缀 ——")
    print(f"  去重后 base 反应 {d['n_base']}（酶促 {d['n_enzymatic']}）："
          f"kcat 覆盖 {d['kcat_coverage_of_enzymatic']:.2%} ｜ EC 覆盖 {d['ec_coverage_of_enzymatic']:.2%}")
    print(f"  ① 完全暗物质 {d['dark_full']} 条（占酶促 {d['dark_frac_of_enzymatic']:.2%}）"
          f" ｜ ③ 有kcat·无EC {d['ec_missing_only']} 条（注释可零成本回填）")
    tb = res["top50_bottleneck"]
    print()
    print(f"  —— Top50 瓶颈酶清单（酶池占用 62.81%）交叉 ——")
    print(f"  无 EC 号条目: {tb['n_without_ec']}/{tb['n_listed']}，合计酶池占用 {tb['pool_pct_without_ec']}%")
    if args.json:
        print(json.dumps(res, ensure_ascii=False, indent=2))
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    json.dump(res, open(OUT_JSON, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print(f"\n→ 证据已写入 {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
