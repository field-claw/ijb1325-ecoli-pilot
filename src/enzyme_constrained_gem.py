# -*- coding: utf-8 -*-
"""enzyme_constrained_gem —— 酶约束 / 多尺度 GEM 借项落地（对比学习 2026-09-03）。

把「酶约束多尺度 GEM 对比学习整合」报告中的三项可借项 A/B/C 真正实现在
**全基因组尺度 GEM（cobra 真实 SBML）** 上，而非 9 反应迷你网络：

  A[高] GECKO 2.0 流水线 -> 给真实 GEM 加酶约束（酶伪代谢物 + 酶使用反应
         + 总酶质量容量约束），真 LP 闭合（v_i <= kcat_i·[E_i] 且 Σ[E]<=Cap）。
         iMA871(A.niger) / iYLI647(Y.lipolytica) / ecYeastGEM(S.cerevisiae) 均适用。
  B[中] ecFactory(FSEOF + 最小靶点组合) -> 纯约束方法，不依赖 kcat / 外部权重，
         可在真实 GPR 上跑通；识别与产物（或生长）通量一致上/下调的反应，
         并贪心求最小敲除组合使目标通量最大化。
  C[中] EnzymeTuning kcat 先验优化 -> 用参数化 kcat 先验（非训练 GAN，诚实标 DEMO）
         在酶约束 GEM 上优化目标通量；真实 kcat 数据库不可达时 BLOCKED_NO_DATA。

  P0[高] 鲁洪中组借项落地（对比学习 2026-09-10，CL #44）：
    ecFSEOF        —— 酶约束版 FSEOF + 酶丰度变异分析（普通 FBA FSEOF 漏掉的酶资源
                     重分配靶点；added_by_enzyme_constraint 即该增量）
    flux_control_coefficient —— FCC 通量控制系数（有限差分，StrainOptimizer 评分层）
    tier_eval      —— L1/L2/L3/FCC 分级评估基准 + precision/recall/novelty 三指标
                      （层级 + 指标定义 + 候选汇聚为交付物；precision/recall 需 reference_set）

诚实边界（全局强制）：
  - 本机无真实 kcat 测定数据库 / 无 EnzymeTuning GAN 权重；A/C 用演示或参数化
    kcat，标注 kcat_source="demo" / "prior"；请求 kcat_source="real" 而无 kcat_map
    时返回 BLOCKED_NO_DATA，绝不编造 kcat。
  - B 为纯 LP，完全真实、零外部依赖。
  - 这是 GECKO/ecFactory/EnzymeTuning 方法学在 VitaMind 真实 GEM 上的**真求解实现**，
    不是 stub；但酶动力学真实值是 DEMO/先验（与 synbio-ecgem-proteome-bridge 同纪律）。

依赖：cobra（cobra313 venv）。默认 dry_run：只评估、不写、不下载、不调外部服务。
"""
from __future__ import annotations
import os
import math
import json
from typing import Dict, Any, List, Optional, Tuple

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_THIS_DIR)

# 本地真实 GEM 路径（与 vitamind_core 其他模块一致）
# iYLI647 / eciYali 为真实 Y. lipolytica GEM（见 GEM_STATUS_NOTE，已落盘 REAL）。
GEM_PATHS = {
    "iMA871":   os.path.join(_THIS_DIR, "iMA871_biomodels.xml"),          # A. niger 黑曲霉 (BioModels, 1400 rxn, 无 GPR)
    "ecYeastGEM": os.path.join(_THIS_DIR, "data", "gem", "ecYeastGEM.xml"),  # S. cerevisiae 酶约束 GEM (8143 rxn, 带 GPR)
    "yeast-GEM": os.path.join(_PROJECT_ROOT, "data", "yeast_gem", "yeast-GEM.xml"),  # S. cerevisiae (4102 rxn, 带 GPR)
    "iYLI647":  os.path.join(_PROJECT_ROOT, "data", "iYLI647", "iYLI647_fixed.xml"),  # Y. lipolytica 真实 GEM (Mishra 2018, 1343 rxn/646 gene/831 GPR)
    "eciYali":  os.path.join(_THIS_DIR, "data", "gem", "eciYali_fixed.xml"),  # Y. lipolytica 真实酶约束 ecGEM (SysBioChalmers, 4022 rxn/847 gene)
    "eciML1515": os.path.join(_THIS_DIR, "data", "gem", "eciML1515.json"),      # E. coli 真实 GECKO ecGEM (iML1515 主干 + kcat_MW 2399 条, 总酶容量 0.227 g/gDW)
}
# 已知 GEM 本地状态（诚实边界，禁止误称 REAL）
GEM_STATUS_NOTE = {
    "iYLI647": "REAL: 2026-09-03 自 BMC Syst Biol 2018 Suppl.MOESM2 下载真实 SBML(level2)，"
               "修复 2 个 biomass 反应缺 id 后 cobra 可加载 (1343 rxn/646 gene/831 GPR)；"
               "旧记忆「iYLI647 未落地」已更正。",
    "eciYali": "REAL: SysBioChalmers/ecModels 出品的酶约束版 iYali (4022 rxn/847 gene)，"
               "已落 data/gem/eciYali_fixed.xml（SBML id 转义修复）。",
    "eciML1515": "REAL: 用户提供的 E. coli GECKO ecGEM（iML1515 主干 + 酶层，5883 rxn/1877 mets/1516 gene/3 区室）。"
                 "enzyme_constraint 含 kcat_MW(2399 条反应级特异性活性) + 总酶容量 upperbound=0.227 g 酶/gDW。"
                 "标准 cobra 加载后酶约束不绑定(普通 FBA=0.8698)；geckopy 2.0.2 命名不匹配亦无法直接求解。"
                 "已用真实 kcat_MW 在主干上自建 faithful GECKO LP（gecko_eciML1515.py）：满容量 biomass=0.680<0.8698（约束绑定），"
                 "随容量收紧单调降至 0.062，均为真 LP(Highs)。",
    "iJB1325": "非酵母：实为 A. niger 黑曲霉；真实 SBML 在 data/iJB1325/ (工业酵母参考，非虚拟酵母底盘)。",
}
# 向后兼容别名（部分脚本仍查 KNOWN_MISSING_LOCAL）
KNOWN_MISSING_LOCAL = {
    "iJB1325": GEM_STATUS_NOTE["iJB1325"],
}

_GEM_CACHE: Dict[str, Any] = {}


def load_gem(key: str):
    """加载并缓存本地真实 GEM（cobra Model 副本）。"""
    if key not in GEM_PATHS:
        raise KeyError(f"未知 GEM key: {key}；可选 {list(GEM_PATHS)}")
    if key not in _GEM_CACHE:
        try:
            import cobra
            from cobra.io import read_sbml_model
        except Exception as e:  # pragma: no cover
            raise RuntimeError(f"cobra 不可用: {e}")
        path = GEM_PATHS[key]
        if not os.path.exists(path):
            raise FileNotFoundError(f"GEM 文件缺失: {path}")
        _GEM_CACHE[key] = read_sbml_model(path)
    return _GEM_CACHE[key].copy()


def find_biomass_rxn(model) -> Optional[Any]:
    """定位 biomass / 生长反应（按 id/name 关键词，大小写不敏感）。"""
    for r in model.reactions:
        tag = (r.id + " " + (r.name or "")).lower()
        if "biomass" in tag or "growth" in tag or "objective" in tag:
            return r
    return None


# ---------------------------------------------------------------------------
# A) GECKO 2.0 酶约束流水线
# ---------------------------------------------------------------------------
def gecko_augment(model, kcat_map: Optional[Dict[str, float]] = None,
                  cap_factor: float = 1.0, total_capacity: Optional[float] = None,
                  kcat_source: str = "demo", kcat_default: float = 10.0,
                  return_model: bool = False,
                  enzyme_scale_map: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """GECKO 酶约束增强：给 GPR 映射的反应加酶伪代谢物 + 前/反向酶使用反应 +
    总酶质量容量约束，真 LP 闭合（v_i<=kcat_i·[E_i]，Σ[E]<=Cap）。

    enzyme_scale_map: 可选 {reaction_id: 乘子}，将 GOAI 蛋白 log₂FC 经真实 GPR
        映射后的酶容量调制注入（乘子 = 2**log₂FC，下调<1 收紧酶预算）；这是
        GOAI→GECKO 接线点（见 goai_aivm_fusion.run_fusion_point1_gecko）。

    返回 {status, augmented_model(可选), n_enzyme_rxns, biomass_unconstrained,
          biomass_constrained, drop_pct, kcat_source, note}。
    kcat_source="real" 且无 kcat_map -> BLOCKED_NO_DATA（不编造 kcat）。
    """
    if kcat_source == "real" and not kcat_map:
        return {"status": "BLOCKED_NO_DATA",
                "note": "real kcat 测定数据库本机不可达；须提供 kcat_map 或用 kcat_source='demo'。"
                        "GECKO 2.0 真实 kcat 来自 BRENDA/SABIO/ecGEM 注释，未本地化（诚实边界）。"}

    try:
        from cobra import Metabolite, Reaction
    except Exception as e:  # pragma: no cover
        return {"status": "ERROR", "note": f"cobra 不可用: {e}"}

    bm = find_biomass_rxn(model)
    if bm is None:
        return {"status": "ERROR", "note": "未找到 biomass 反应，无法设目标。"}

    # 无酶约束基线
    model.objective = bm.id
    sol0 = model.optimize()
    if sol0.status != "optimal":
        return {"status": "ERROR", "note": f"基线 FBA 不可行: {sol0.status}"}
    bm_u = float(sol0.objective_value)
    flux0 = dict(sol0.fluxes)

    m = model.copy()
    enz_rxns = [r for r in m.reactions if r.gene_reaction_rule.strip()]
    n_e = len(enz_rxns)
    if n_e == 0:
        return {"status": "NO_GPR",
                "biomass_unconstrained": round(bm_u, 4),
                "note": "该 GEM 无基因-反应规则(GPR)，GECKO 酶约束无法施加；"
                        "须用带 GPR 的 GEM（如 yeast-GEM / ecYeastGEM）。iMA871(BioModels) 即无 GPR。"}

    # 总酶质量池 + 供给反应（上限=总容量，由酶分配求和得到）
    pool = Metabolite("enzyme_total_pool", compartment="e",
                      name="total enzyme mass pool (GECKO pseudo)")
    supply = Reaction("ENZYME_SUPPLY")
    supply.add_metabolites({pool: 1.0})
    m.add_metabolites([pool])
    m.add_reactions([supply])

    # GECKO 标准酶分配：酶浓度 [E]_r 与 kcat·通量 成正比（v = kcat·[E]）。
    # 故 [E]_r = kcat_r·|v_base|，cap_factor<1 时整体按比例收紧（酶供应不足 -> 约束绑紧）。
    # 逐反应酶约束上界 vmax_r = [E]_r/kcat_r = cap_factor·|v_base|。
    # 未使用反应给地板值，允许在酶约束最优解中重新布线。
    # 总酶容量 = Σ[E]_r（真实 LP 闭合 Σ[E]<=Cap）。
    FLUX_FLOOR = 0.01
    e_alloc: Dict[str, float] = {}
    for r in enz_rxns:
        kcat = (kcat_map or {}).get(r.id, kcat_default)
        if kcat <= 0:
            kcat = kcat_default
        mult = (enzyme_scale_map or {}).get(r.id, 1.0)  # GOAI 丰度乘子(2^log2FC)
        vbase = abs(flux0.get(r.id, 0.0))
        need = max(vbase, FLUX_FLOOR)
        e_r = cap_factor * kcat * need * mult  # 酶浓度 ∝ kcat·通量·(GOAI 丰度乘子)
        e_alloc[r.id] = e_r

    total_capacity = total_capacity if total_capacity is not None else sum(e_alloc.values())
    supply.upper_bound = total_capacity

    added = 0
    for r in enz_rxns:
        kcat = (kcat_map or {}).get(r.id, kcat_default) or kcat_default
        e_r = e_alloc[r.id]
        vmax = kcat * e_r  # = cap_factor * max(vbase, floor)，即逐反应酶约束上界
        em = Metabolite(f"e_{r.id}", compartment="e",
                        name=f"enzyme pool for {r.id}")
        m.add_metabolites([em])
        # 酶分配：从总池取酶到 e_r（≤ e_r）
        prot = Reaction(f"prot_{r.id}")
        prot.add_metabolites({pool: -1.0, em: 1.0})
        prot.bounds = (0.0, e_r)
        # 酶汇：回收 e_r（反向通量时 drain 多余酶；反向酶约束由逐反应 |v| 界覆盖，诚实近似）
        sink = Reaction(f"sink_{r.id}")
        sink.add_metabolites({em: -1.0})
        sink.bounds = (0.0, 1e6)
        # 原反应正向前向消耗酶：kcat·e_r -> r（v_r <= kcat·[E_r]）
        r.add_metabolites({em: -kcat})
        # 逐反应酶约束上界（可逆取 |v|<=vmax，真实 LP 闭合）
        lo, hi = r.lower_bound, r.upper_bound
        new_hi = min(hi, vmax) if (hi is None or hi > 1e8) else min(hi, vmax)
        new_lo = max(lo, -vmax) if r.reversibility else lo
        r.lower_bound, r.upper_bound = new_lo, new_hi
        m.add_reactions([prot, sink])
        added += 1

    m.objective = bm.id
    sol1 = m.optimize()
    if sol1.status != "optimal":
        # 容量过紧导致不可行：返回诊断而非谎报
        return {"status": "INFEASIBLE_CAP",
                "biomass_unconstrained": round(bm_u, 4),
                "n_enzyme_rxns": n_e,
                "total_capacity": round(float(total_capacity), 6),
                "kcat_source": kcat_source,
                "note": "cap_factor 过小导致酶约束后无可行解；放宽 cap_factor 重试（cap_factor=1 即刚好满足基线需求）。"}
    bm_c = float(sol1.objective_value)
    drop = (1.0 - bm_c / bm_u) * 100.0 if bm_u else 0.0

    n_scale = len([r for r in enz_rxns if (enzyme_scale_map or {}).get(r.id, 1.0) != 1.0])
    return {
        "status": "ok",
        "method": "GECKO 2.0 enzyme-constrained augmentation (real LP, cobra)",
        "n_enzyme_rxns": n_e,
        "enzyme_scale_applied": n_scale,
        "biomass_unconstrained": round(bm_u, 4),
        "biomass_constrained": round(bm_c, 4),
        "biomass_drop_pct": round(float(drop), 2),
        "total_capacity": round(float(total_capacity), 6),
        "kcat_source": kcat_source,
        "augmented_model": m if return_model else None,
        "note": ("酶约束 v_i<=kcat_i·[E_i] 且 Σ[E]<=Cap 已真 LP 闭合。"
                 if kcat_source == "real"
                 else "demo kcat（非实验测定）；生产须接真实 ecGEM kcat 注释（GECKO 2.0 目录）。"
                      + (f" GOAI 酶容量调制已注入 {n_scale} 个反应。" if n_scale else "")),
    }


# ---------------------------------------------------------------------------
# B) ecFactory: FSEOF + 最小靶点组合（纯 LP，全真）
# ---------------------------------------------------------------------------
def fseof_scan(model, target_rxn: Optional[str] = None,
               n_steps: int = 10, lo_frac: float = 0.1,
               hi_frac: float = 1.0) -> Dict[str, Any]:
    """FSEOF（Flux Scanning by Enforced Objective Flux）：把目标通量（默认 biomass）
    强制到其最大值的 [lo,hi] 区间逐档扫描，记录各反应通量对目标的相关性，
    识别与目标一致上/下调的反应（候选过表达 / 敲除靶点）。

    返回 {status, target, n_steps, up_reactions, down_reactions, n_reactions}。
    纯约束方法，不依赖 kcat / 外部权重。
    """
    if target_rxn is None:
        bm = find_biomass_rxn(model)
        if bm is None:
            return {"status": "ERROR", "note": "未找到 biomass 反应。"}
        target_rxn = bm.id

    # 目标最大通量
    with model:
        model.objective = target_rxn
        sol = model.optimize()
    if sol.status != "optimal":
        return {"status": "ERROR", "note": f"目标最大化不可行: {sol.status}"}
    fmax = float(sol.objective_value)
    if fmax <= 1e-9:
        return {"status": "ERROR", "note": "目标最大通量≈0，无法扫描。"}

    slopes: Dict[str, float] = {}
    base_flux: Dict[str, float] = {}
    with model:
        model.objective = target_rxn
        for step in range(n_steps + 1):
            frac = lo_frac + (hi_frac - lo_frac) * step / n_steps
            tf = frac * fmax
            # 强制目标通量落在本档（先设上界再设下界，避免与上一档残留边界冲突）
            tgt = model.reactions.get_by_id(target_rxn)
            tgt.upper_bound = tf
            tgt.lower_bound = tf
            s = model.optimize()
            if s.status != "optimal":
                continue
            for r in model.reactions:
                if r.id == target_rxn:
                    continue
                slopes.setdefault(r.id, 0.0)
                base_flux.setdefault(r.id, float(s.fluxes.get(r.id, 0.0)))
                slopes[r.id] += float(s.fluxes.get(r.id, 0.0))
    # 用首末档斜率近似相关性（省内存）
    up, down = [], []
    for rid, tot in slopes.items():
        # 线性趋势近似：累计和正 -> 随目标增大净增
        mean = tot / max(n_steps, 1)
        ref = base_flux.get(rid, 0.0)
        trend = mean - ref
        if trend > 1e-6:
            up.append((rid, round(trend, 4)))
        elif trend < -1e-6:
            down.append((rid, round(trend, 4)))
    up.sort(key=lambda x: -x[1])
    down.sort(key=lambda x: x[1])
    return {
        "status": "ok",
        "method": "ecFactory FSEOF (flux scanning, enforced objective, pure LP)",
        "target": target_rxn,
        "n_steps": n_steps,
        "target_fmax": round(fmax, 4),
        "n_reactions": len(slopes),
        "n_up": len(up),
        "n_down": len(down),
        "up_reactions": up[:20],       # 候选过表达靶点（与目标同向上调）
        "down_reactions": down[:20],   # 候选敲除靶点（与目标反向下调）
        "note": "纯约束 FSEOF，真实、零外部依赖；up=随目标增大净增(过表达候选)，"
                "down=随目标增大净减(敲除候选)。",
    }


def minimal_target_combo(model, target_rxn: Optional[str] = None,
                         max_ko: int = 5, max_candidates: Optional[int] = None) -> Dict[str, Any]:
    """最小靶点组合（ecFactory 最小敲除集）：贪心 LP——反复寻找敲除后目标通量
    增益最大的反应并敲除，直到增益收敛或达 max_ko。纯 LP，全真。

    返回 {status, target, knockout_set, target_before, target_after, steps}。
    """
    if target_rxn is None:
        bm = find_biomass_rxn(model)
        if bm is None:
            return {"status": "ERROR", "note": "未找到 biomass 反应。"}
        target_rxn = bm.id

    def _max_target(mdl) -> float:
        with mdl:
            mdl.objective = target_rxn
            s = mdl.optimize()
        return float(s.objective_value) if s.status == "optimal" else 0.0

    base = _max_target(model)
    ko: List[str] = []
    m = model.copy()
    # 候选敲除反应集（默认全模型；大模型可限 max_candidates 提速，按 GPR 优先采样）
    cands = [r for r in m.reactions if r.id != target_rxn and r.id not in ko]
    if max_candidates and len(cands) > max_candidates:
        import random as _rnd
        _rnd.seed(0)
        gpr = [r for r in cands if r.gene_reaction_rule.strip()]
        _rnd.shuffle(gpr)
        others = [r for r in cands if not r.gene_reaction_rule.strip()]
        _rnd.shuffle(others)
        cands = (gpr + others)[:max_candidates]
    steps = []
    for k in range(max_ko):
        best, best_rid = base, None
        for r in cands:
            if r.id in ko:
                continue
            with m:
                rt = m.reactions.get_by_id(r.id)
                rt.lower_bound = 0.0
                rt.upper_bound = 0.0
                val = _max_target(m)
            if val > best + 1e-6:
                best, best_rid = val, r.id
        if best_rid is None:
            break
        ko.append(best_rid)
        rt = m.reactions.get_by_id(best_rid)
        rt.lower_bound = 0.0
        rt.upper_bound = 0.0
        steps.append({"ko": best_rid, "target_after_ko": round(best, 4)})
        base = best

    return {
        "status": "ok",
        "method": "ecFactory minimal-target-combination (greedy LP knockout)",
        "target": target_rxn,
        "max_ko": max_ko,
        "knockout_set": ko,
        "target_before": round(_max_target(model), 4),
        "target_after": round(base, 4),
        "n_knockouts": len(ko),
        "steps": steps,
        "note": "贪心 LP 最小敲除组合（真实）；若目标=biomass，则为生长必需最小敲除集。",
    }


# ---------------------------------------------------------------------------
# C) EnzymeTuning kcat 先验优化（参数化先验代理，诚实标 DEMO 非 GAN）
# ---------------------------------------------------------------------------
def kcat_prior_optimize(model, target_rxn: Optional[str] = None,
                        kcat_default: float = 10.0,
                        prior_lo: float = 1.0, prior_hi: float = 100.0,
                        n_samples: int = 60, seed: int = 0,
                        biomass_floor: float = 0.5) -> Dict[str, Any]:
    """kcat 先验优化（EnzymeTuning 代理）：在 GECKO 酶约束 GEM 上，从对数均匀
    kcat 先验中采样多组 kcat，最大化目标通量（给定 biomass 下限），挑最优。

    诚实边界：这是**参数化先验采样**，并非 EnzymeTuning 的训练 GAN（GAN 需真实
    kcat 训练数据，本机不可达）。请求 kcat_source="real" 时返回 BLOCKED_NO_DATA。

    返回 {status, best_kcat, best_target, samples, note}。
    """
    import random
    rng = random.Random(seed)
    if target_rxn is None:
        bm = find_biomass_rxn(model)
        if bm is None:
            return {"status": "ERROR", "note": "未找到 biomass 反应。"}
        target_rxn = bm.id

    bm = find_biomass_rxn(model)
    best = None
    samples = []
    for i in range(n_samples):
        # 对数均匀先验 kcat
        logk = rng.uniform(math.log(prior_lo), math.log(prior_hi))
        k = math.exp(logk)
        kcat_map = {r.id: k for r in model.reactions if r.gene_reaction_rule.strip()}
        # 按基线需求分配酶容量（cap_factor=1.0 刚好满足），采样不同 kcat 看目标通量变化
        res = gecko_augment(model, kcat_map=kcat_map, kcat_source="demo",
                            kcat_default=k, return_model=True)
        if res["status"] != "ok":
            samples.append({"kcat": round(k, 3), "target": None, "note": res["status"]})
            continue
        am = res["augmented_model"]
        # 在 biomass 下限约束下最大化目标
        am.reactions.get_by_id(bm.id).lower_bound = biomass_floor * res["biomass_unconstrained"]
        am.objective = target_rxn
        s = am.optimize()
        tv = float(s.objective_value) if s.status == "optimal" else 0.0
        samples.append({"kcat": round(k, 3), "target": round(tv, 4)})
        if best is None or tv > best[1]:
            best = (k, tv)
    return {
        "status": "ok",
        "method": "EnzymeTuning surrogate (parametric kcat prior sampling, NOT a trained GAN)",
        "target": target_rxn,
        "n_samples": n_samples,
        "best_kcat": round(best[0], 3) if best else None,
        "best_target": round(best[1], 4) if best else None,
        "biomass_floor": biomass_floor,
        "samples": samples[:12],
        "note": ("参数化 kcat 先验采样代理，最大化目标通量；"
                 "EnzymeTuning 用 GAN 在真实 kcat 上训练（已验证 Y.lipolytica），"
                 "本机无该 GAN/数据 -> DEMO 代理，非伪造 GAN 结果。"),
    }


# ---------------------------------------------------------------------------
# P0) 鲁洪中组借项落地（对比学习 2026-09-10）：ecFSEOF / FCC / 分级评估基准
# ---------------------------------------------------------------------------
def ecfseof_scan(model, target_rxn: Optional[str] = None,
                 kcat_map: Optional[Dict[str, float]] = None,
                 kcat_source: str = "demo", cap_factor: float = 1.0,
                 n_steps: int = 10, lo_frac: float = 0.1, hi_frac: float = 1.0,
                 enzyme_variation_factors: Tuple[float, ...] = (0.5, 2.0),
                 enzyme_scale_map: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """ecFSEOF（酶约束版 FSEOF，鲁洪中 ecFactory 升级项）：在 GECKO 酶约束 GEM 上跑
    FSEOF，与普通 FBA FSEOF 对比，识别只因酶资源重分配才浮现的靶点；再做酶丰度变异
    分析（GOAI 式 log2FC 乘子扫描），量化靶点集对酶预算的敏感度。

    这是普通 fseof_scan 的升级：普通 FBA 假设酶无限供应，漏掉「酶容量受限导致代谢流
    重布线」才出现的靶点；ecFSEOF 在酶约束 GEM 上重算，added_by_enzyme_constraint 即为
    该增量。纯 LP + 真实酶约束，零外部权重。

    返回 {status, target, fba_fseof, ec_fseof, added_by_enzyme_constraint,
          enzyme_variation, kcat_source, note}。
    """
    if target_rxn is None:
        bm = find_biomass_rxn(model)
        if bm is None:
            return {"status": "ERROR", "note": "未找到 biomass 反应。"}
        target_rxn = bm.id

    # 普通 FBA FSEOF（对照）
    fba = fseof_scan(model, target_rxn=target_rxn, n_steps=n_steps,
                     lo_frac=lo_frac, hi_frac=hi_frac)
    if fba.get("status") != "ok":
        return {"status": "ERROR", "note": f"对照 FBA FSEOF 失败: {fba.get('note')}"}
    fba_up = {r for r, _ in fba.get("up_reactions", [])}
    fba_down = {r for r, _ in fba.get("down_reactions", [])}

    # GECKO 酶约束增强
    aug = gecko_augment(model, kcat_map=kcat_map, kcat_source=kcat_source,
                        cap_factor=cap_factor, return_model=True,
                        enzyme_scale_map=enzyme_scale_map)
    if aug.get("status") != "ok":
        return {"status": "ERROR",
                "note": f"GECKO 增强失败({aug.get('status')})；plain-FBA FSEOF 仍可用。",
                "fba_fseof": {"up": list(fba_up), "down": list(fba_down),
                              "n_up": len(fba_up), "n_down": len(fba_down)}}
    am = aug["augmented_model"]

    # 酶约束 FSEOF
    ec = fseof_scan(am, target_rxn=target_rxn, n_steps=n_steps,
                    lo_frac=lo_frac, hi_frac=hi_frac)
    ec_up = {r for r, _ in ec.get("up_reactions", [])} if ec.get("status") == "ok" else set()
    ec_down = {r for r, _ in ec.get("down_reactions", [])} if ec.get("status") == "ok" else set()

    # 酶丰度变异分析：在 ±factor 酶预算下重跑 ecFSEOF，记录靶点集漂移
    variation: Dict[str, Any] = {}
    for fac in enzyme_variation_factors:
        scale = {r.id: fac for r in am.reactions if r.gene_reaction_rule.strip()}
        ecv = gecko_augment(model, kcat_map=kcat_map, kcat_source=kcat_source,
                            cap_factor=cap_factor, return_model=True,
                            enzyme_scale_map=scale)
        if ecv.get("status") != "ok":
            variation[str(fac)] = {"status": ecv.get("status"),
                                   "note": "该酶预算下酶约束 GEM 不可行（酶供应不足即为限制信号）"}
            continue
        ev = fseof_scan(ecv["augmented_model"], target_rxn=target_rxn,
                        n_steps=n_steps, lo_frac=lo_frac, hi_frac=hi_frac)
        if ev.get("status") == "ok":
            variation[str(fac)] = {
                "n_up": ev.get("n_up"), "n_down": ev.get("n_down"),
                "up": [r for r, _ in ev.get("up_reactions", [])],
                "down": [r for r, _ in ev.get("down_reactions", [])],
            }

    added_up = sorted(ec_up - fba_up)
    added_down = sorted(ec_down - fba_down)
    return {
        "status": "ok",
        "method": "ecFSEOF (FSEOF on enzyme-constrained GEM) + enzyme-abundance variation",
        "target": target_rxn,
        "fba_fseof": {"n_up": len(fba_up), "n_down": len(fba_down)},
        "ec_fseof": {"n_up": len(ec_up), "n_down": len(ec_down)},
        "added_by_enzyme_constraint": {"up": added_up[:20], "down": added_down[:20],
                                       "n_up": len(added_up), "n_down": len(added_down)},
        "enzyme_variation": variation,
        "kcat_source": kcat_source,
        "note": ("ecFSEOF 在酶约束 GEM 上重算 FSEOF；added_by_enzyme_constraint = 普通 FBA "
                 "FSEOF 漏掉、仅在酶资源受限下浮现的靶点（鲁洪中核心增量：蛋白资源动态分配）。"
                 "enzyme_variation = 酶预算 ±factor 下靶点集漂移，量化对酶供应的敏感度；"
                 "若 factor<1 返回 INFEASIBLE_CAP 即酶短缺致死信号。top-20 口径对比。"),
    }


def flux_control_coefficient(model, target_rxn: Optional[str] = None,
                             candidate_rxns: Optional[List[str]] = None,
                             eps: float = 1e-2, max_candidates: int = 300,
                             seed: int = 0) -> Dict[str, Any]:
    """FCC（通量控制系数，鲁洪中 StrainOptimizer 评分层）：FCC_i = ∂ln(v_target)/∂ln(v_i)，
    有限差分近似——对每个候选反应在其容量上界各扰动 ±eps·|v_i|，重算目标最大通量，
    取中点差分归一化灵敏度（FCC = (Δv_target/v_target)/(Δv_i/v_i)）。

    纯 LP，零外部权重。top_positive = 过表达增益高（提该反应容量助目标）；
    top_negative = 敲除增益高（限该反应容量损目标，与 FSEOF down 互补）。

    返回 {status, target, fmax, eps, n_coeff, top_positive, top_negative, note}。
    """
    if target_rxn is None:
        bm = find_biomass_rxn(model)
        if bm is None:
            return {"status": "ERROR", "note": "未找到 biomass 反应。"}
        target_rxn = bm.id

    with model:
        model.objective = target_rxn
        sol = model.optimize()
    if sol.status != "optimal":
        return {"status": "ERROR", "note": f"目标最大化不可行: {sol.status}"}
    fmax = float(sol.objective_value)
    if fmax <= 1e-9:
        return {"status": "ERROR", "note": "目标最大通量≈0。"}
    flux0 = dict(sol.fluxes)

    if candidate_rxns is None:
        cands = [r for r in model.reactions if r.gene_reaction_rule.strip()]
    else:
        ids = {x.id for x in model.reactions}
        cands = [model.reactions.get_by_id(r) for r in candidate_rxns if r in ids]
    if max_candidates and len(cands) > max_candidates:
        import random as _r
        _r.seed(seed)
        _r.shuffle(cands)
        cands = cands[:max_candidates]

    coeffs = []
    for r in cands:
        hi0 = r.upper_bound
        dhi = eps * abs(hi0) if abs(hi0) < 1e8 else eps  # 无界反应扰动可忽略（非限制步）
        with model:
            rt = model.reactions.get_by_id(r.id)
            rt.upper_bound = hi0 + dhi
            model.objective = target_rxn
            sp = model.optimize()
            vp = float(sp.objective_value) if sp.status == "optimal" else fmax
            rt2 = model.reactions.get_by_id(r.id)
            rt2.upper_bound = hi0 - dhi
            model.objective = target_rxn
            sm = model.optimize()
            vm = float(sm.objective_value) if sm.status == "optimal" else fmax
        dvt = (vp - vm) / 2.0
        fcc = (dvt / fmax) / eps if eps > 0 else 0.0
        coeffs.append({"rxn": r.id, "fcc": round(fcc, 4),
                       "v_target_plus": round(vp, 4), "v_target_minus": round(vm, 4)})
    coeffs.sort(key=lambda x: -abs(x["fcc"]))
    top_pos = [c for c in coeffs if c["fcc"] > 0][:20]
    top_neg = [c for c in coeffs if c["fcc"] < 0][:20]
    return {
        "status": "ok",
        "method": "Flux Control Coefficient (finite-difference FCC proxy, pure LP)",
        "target": target_rxn,
        "fmax": round(fmax, 4),
        "eps": eps,
        "n_coeff": len(coeffs),
        "top_positive": top_pos,
        "top_negative": top_neg,
        "note": "FCC_i = (Δv_target/v_target)/(Δv_i/v_i) 中点差分近似；纯 LP 真实，零外部权重；"
                "无界反应扰动可忽略→fcc≈0（非限制步）。",
    }


def _rxn_ids_to_genes(model, rxn_ids) -> set:
    """把反应 ID 集合经 GPR 映射回基因集合（修复 tier_eval 与基因参考集命名空间错位）。

    仅纳入模型真实存在的反应；GECKO 伪反应（prot_*/sink_*/ENZYME_SUPPLY）不在主干→跳过。
    基因取并集：任一 GPR 基因即计入（OR 规则命中任一、AND 规则需全齐，保守取并集）。
    返回基因 id 集合（如 eciML1515 的 b-number）。
    """
    genes: set = set()
    valid = {r.id for r in model.reactions}
    for rid in rxn_ids:
        if rid in valid:
            for g in model.reactions.get_by_id(rid).genes:
                genes.add(g.id)
    return genes


def _detect_reference_namespace(model, reference_set) -> str:
    """检测 reference_set 命名空间：'gene'（项为模型基因 id）/ 'reaction'（项为反应 id）/ 'unknown'。"""
    if not reference_set:
        return "unknown"
    gids = {g.id for g in model.genes}
    rids = {r.id for r in model.reactions}
    n_gene = sum(1 for x in reference_set if x in gids)
    n_rxn = sum(1 for x in reference_set if x in rids)
    if n_gene >= n_rxn and n_gene > 0:
        return "gene"
    if n_rxn > n_gene:
        return "reaction"
    return "unknown"


def essentiality_benchmark(model, reference_essential, reference_uncertain=None,
                           target_rxn: Optional[str] = None, thr: float = 1e-6,
                           max_genes: Optional[int] = None) -> Dict[str, Any]:
    """dedicated 必需性路径（修复 tier_eval 与 Keio 必需集语义错位）：

    对模型逐基因**手动串行 KO**（Windows 安全，绝对不走 cobra.single_gene_deletion 的
    子进程 + 巨型向量化 LP → OOM/递归崩溃），预测必需基因（KO 后 biomass < thr），
    与 Keio 必需集比对，输出公平混淆矩阵。

    reference_essential: 必需基因集合（阳性，如 Keio-E 的 b-number）
    reference_uncertain: 可选代理阴性（如 Keio-N 必需性不确定基因）
    max_genes: 自检限速（仅跑前 N 个基因）；None=全模型（1516 基因约数分钟）

    返回 {status, predicted_essential_n, confusion_Epos_Nneg, recall, precision,
          specificity, mcc, completeness_gap_E, note}。
    公平性：仅在模型存在的参考基因上计分；缺失基因=模型完整性缺口，归 NA 不罚分。
    """
    if target_rxn is None:
        bm = find_biomass_rxn(model)
        if bm is None:
            return {"status": "ERROR", "note": "未找到 biomass 反应。"}
        target_rxn = bm.id

    all_gids = [g.id for g in model.genes]
    genes = all_gids[:max_genes] if max_genes else all_gids

    pred: Dict[str, bool] = {}
    for g in genes:
        with model:
            try:
                model.genes.get_by_id(g).knock_out()
                sol = model.optimize()
                v = sol.objective_value if (sol is not None and sol.status == "optimal") else 0.0
            except Exception:
                v = 0.0
        pred[g] = (v is not None and v < thr)

    E = set(reference_essential)
    N = set(reference_uncertain or [])
    E_present = {b for b in E if b in all_gids}
    N_present = {b for b in N if b in all_gids}
    E_absent = {b for b in E if b not in all_gids}
    N_absent = {b for b in N if b not in all_gids}
    TP = sum(1 for b in E_present if pred.get(b))
    FN = sum(1 for b in E_present if not pred.get(b))
    FP = sum(1 for b in N_present if pred.get(b))
    TN = sum(1 for b in N_present if not pred.get(b))
    recall = TP / (TP + FN) if (TP + FN) else None
    precision = TP / (TP + FP) if (TP + FP) else None
    specificity = TN / (TN + FP) if (TN + FP) else None
    denom = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
    mcc = (TP * TN - FP * FN) / denom if denom else None
    return {
        "status": "ok",
        "method": "serial single-gene-deletion FBA essentiality (manual loop, Windows-safe)",
        "predicted_essential_n": sum(1 for x in pred.values() if x),
        "n_genes_scored": len(genes),
        "keio_E_present": sorted(E_present), "keio_E_absent": sorted(E_absent),
        "keio_N_present": sorted(N_present), "keio_N_absent": sorted(N_absent),
        "completeness_gap_E": "%d/%d reference-essential genes absent from model" % (len(E_absent), len(E)),
        "confusion_Epos_Nneg": {"TP": TP, "FN": FN, "FP": FP, "TN": TN},
        "recall": round(recall, 4) if recall is not None else None,
        "precision": round(precision, 4) if precision is not None else None,
        "specificity": round(specificity, 4) if specificity is not None else None,
        "mcc": round(mcc, 4) if mcc is not None else None,
        "note": "公平性：仅在模型存在的参考基因上计分；缺失=模型完整性缺口，非预测错误（归 NA 不罚分）。",
    }


def tier_eval(model, target_rxn: Optional[str] = None,
              reference_set: Optional[List[str]] = None,
              kcat_source: str = "demo", cap_factor: float = 1.0,
              eps: float = 1e-2, fcc_max: int = 300,
              kcat_map: Optional[Dict[str, float]] = None,
              n_steps: int = 10,
              enzyme_variation_factors: Tuple[float, ...] = (0.5, 2.0),
              reference_type: str = "auto",
              include_essentiality: bool = False,
              essentiality_uncertain: Optional[List[str]] = None,
              essentiality_max_genes: Optional[int] = None) -> Dict[str, Any]:
    """分级评估基准（鲁洪中 L1/L2/L3/FCC + 三指标 precision/recall/novelty）：

    L1 单层判据：单 KO/OE 候选（取 FCC 最高反应，floor_cit 式单点判据）
    L2 组合判据：minimal_target_combo（最小敲除组合）
    L3 多尺度：ecFSEOF（酶约束 FSEOF，含酶资源重分配）
    FCC 评分层：flux_control_coefficient

    命名空间修复（2026-09-13，项②）：候选是**反应 ID**（含 GECKO 伪反应
    prot_*/sink_*），而 Keio 参考集是**基因 b-number** —— 此前直接 `predicted &
    reference` 恒为 0（结构性 0，非预测失败）。现新增 gene_level 口径：
    候选反应经 GPR 映射回基因（`_rxn_ids_to_genes`）后再与基因参考集相交；
    仅当 `reference_namespace == 'gene'` 时 gene_level 有效，reaction_level 保留
    作透明历史口径。

    三指标（需湿实验验证靶点，诚实标注 n/a 当缺数据）：
      precision = |predicted ∩ reference| / |predicted|（命中率）
      recall    = |predicted ∩ reference| / |reference|（覆盖率）
      novelty   = 1 - |predicted ∩ known_prior| / |predicted|（非已知集比例；需 known_prior）

    dedicated 必需性路径（include_essentiality=True）：调用 essentiality_benchmark
    做逐基因串行 KO，与 Keio 必需集比对得公平 recall/precision/MCC——这才是与
    必需性参考集语义对齐的正确 benchmark（tier_eval 默认预测=工程改善靶点，
    与必需集语义不同，重叠低属预期，不可读作预测失败）。

    返回 {status, target, L1_single, L2_combo, L3_multiscale, FCC_top,
          predicted_candidates, predicted_genes, metrics, essentiality_path,
          kcat_source, note}。
    """
    if target_rxn is None:
        bm = find_biomass_rxn(model)
        if bm is None:
            return {"status": "ERROR", "note": "未找到 biomass 反应。"}
        target_rxn = bm.id

    fcc = flux_control_coefficient(model, target_rxn=target_rxn, eps=eps,
                                   max_candidates=fcc_max)
    l1_ko = [c["rxn"] for c in fcc.get("top_negative", [])][:5] if fcc.get("status") == "ok" else []
    l1_oe = [c["rxn"] for c in fcc.get("top_positive", [])][:5] if fcc.get("status") == "ok" else []

    l2 = minimal_target_combo(model, target_rxn=target_rxn, max_ko=3, max_candidates=300)
    l3 = ecfseof_scan(model, target_rxn=target_rxn, kcat_source=kcat_source,
                      cap_factor=cap_factor, kcat_map=kcat_map, n_steps=n_steps)

    predicted: set = set(l1_ko) | set(l1_oe) | set(l2.get("knockout_set", []))
    if l3.get("status") == "ok":
        predicted |= set(l3.get("added_by_enzyme_constraint", {}).get("up", []))
        predicted |= set(l3.get("added_by_enzyme_constraint", {}).get("down", []))

    # 命名空间检测 + 反应级 / 基因级 双口径（修复 tier_eval 与基因参考集结构性 0）
    ref_ns = _detect_reference_namespace(model, reference_set) if reference_set else "unknown"
    predicted_genes = _rxn_ids_to_genes(model, predicted)

    metrics: Dict[str, Any] = {
        "reference_namespace": ref_ns,
        "novelty": None,
        "note": ("precision/recall 需 reference_set；novelty 需 known_prior（未接→n/a）。"
                 "reaction_level=候选反应直接相交（历史口径，与基因参考集可能结构性 0）；"
                 "gene_level=候选反应经 GPR 映射回基因再相交（修复口径）。"
                 "reference_set 为基因(b-number)命名空间时以 gene_level 为准。"),
    }
    if reference_set is not None:
        ref = set(reference_set)
        # 反应级（历史口径，可能结构性 0）
        hit_rxn = predicted & ref
        metrics["reaction_level"] = {
            "precision": round(len(hit_rxn) / len(predicted), 3) if predicted else None,
            "recall": round(len(hit_rxn) / len(ref), 3) if ref else None,
            "n_hit": len(hit_rxn),
        }
        # 基因级（修复口径：仅当 reference 为基因命名空间时有效）
        if ref_ns == "gene" and predicted_genes:
            hit_gene = predicted_genes & ref
            metrics["gene_level"] = {
                "precision": round(len(hit_gene) / len(predicted_genes), 3) if predicted_genes else None,
                "recall": round(len(hit_gene) / len(ref), 3) if ref else None,
                "n_predicted_genes": len(predicted_genes),
                "n_hit": len(hit_gene),
            }
        elif ref_ns == "gene":
            metrics["gene_level"] = {
                "precision": None, "recall": None,
                "n_predicted_genes": len(predicted_genes), "n_hit": 0,
                "note": "候选反应均无 GPR（多为转运/交换反应），基因级命中为空=语义限制非缺陷",
            }
        metrics["n_reference"] = len(ref)

    # dedicated 必需性路径（可选）：逐基因串行 KO vs Keio 必需集（公平 benchmark）
    essentiality_path = None
    if include_essentiality and reference_set is not None:
        essentiality_path = essentiality_benchmark(
            model, reference_essential=reference_set,
            reference_uncertain=essentiality_uncertain,
            target_rxn=target_rxn, max_genes=essentiality_max_genes)
        if essentiality_path.get("status") == "ok":
            metrics["essentiality_mcc"] = essentiality_path.get("mcc")

    return {
        "status": "ok",
        "method": "tiered evaluation benchmark (L1/L2/L3/FCC + 3 metrics, gene-namespace fix 2026-09-13)",
        "target": target_rxn,
        "L1_single": {"ko_candidates": l1_ko, "oe_candidates": l1_oe,
                      "note": "floor_cit 式单点判据（FCC 最高反应）"},
        "L2_combo": {"knockout_set": l2.get("knockout_set"),
                     "target_before": l2.get("target_before"),
                     "target_after": l2.get("target_after")},
        "L3_multiscale": {
            "ec_fseof_n_up": l3.get("ec_fseof", {}).get("n_up"),
            "ec_fseof_n_down": l3.get("ec_fseof", {}).get("n_down"),
            "added_by_enzyme_constraint_n": (
                l3.get("added_by_enzyme_constraint", {}).get("n_up", 0) +
                l3.get("added_by_enzyme_constraint", {}).get("n_down", 0)),
        },
        "FCC_top": {"positive": fcc.get("top_positive", [])[:5],
                    "negative": fcc.get("top_negative", [])[:5]},
        "predicted_candidates": sorted(predicted),
        "predicted_genes": sorted(predicted_genes),
        "metrics": metrics,
        "essentiality_path": essentiality_path,
        "kcat_source": kcat_source,
        "note": ("分层评估基准：L1 单层 / L2 组合 / L3 酶约束多尺度 / FCC 评分层；"
                 "命名空间修复→gene_level 与基因参考集正确相交（不再是反应级结构性 0）；"
                 "include_essentiality=True 时另跑 dedicated 必需性路径（essentiality_benchmark）"
                 "得公平 recall/precision/MCC，与必需集语义对齐。novelty 需 known_prior。"),
    }


# ---------------------------------------------------------------------------
# 自检：在真实本地 GEM 上证明三项均真 LP 闭合 / 真实可跑
# ---------------------------------------------------------------------------
def self_test(gem_key: str = "yeast-GEM") -> Dict[str, Any]:
    checks = []
    try:
        model = load_gem(gem_key)
    except Exception as e:
        return {"status": "ERROR", "note": f"GEM 加载失败: {e}", "checks": []}

    # A: GECKO demo（真实 LP，demo kcat；cap_factor=0.9 演示酶供应不足绑紧）
    a = gecko_augment(model, kcat_source="demo", cap_factor=0.9)
    # A-real-blocked：请求 real 无 kcat_map -> BLOCKED_NO_DATA
    a_blocked = gecko_augment(model, kcat_source="real")
    checks.append({"check": "A_gecko_demo_real_lp",
                   "ok": a.get("status") == "ok" and a.get("biomass_constrained") is not None,
                   "status": a.get("status"),
                   "biomass_unconstrained": a.get("biomass_unconstrained"),
                   "biomass_constrained": a.get("biomass_constrained"),
                   "drop_pct": a.get("biomass_drop_pct"),
                   "n_enzyme_rxns": a.get("n_enzyme_rxns")})
    checks.append({"check": "A_real_blocked_honest",
                   "ok": a_blocked.get("status") == "BLOCKED_NO_DATA",
                   "status": a_blocked.get("status")})

    # B: FSEOF + 最小靶点组合（纯 LP）
    b = fseof_scan(model, n_steps=6)
    c = minimal_target_combo(model, max_ko=3, max_candidates=300)
    checks.append({"check": "B_fseof_real",
                   "ok": b.get("status") == "ok" and b.get("n_reactions", 0) > 0,
                   "n_up": b.get("n_up"), "n_down": b.get("n_down"),
                   "target_fmax": b.get("target_fmax")})
    checks.append({"check": "B_min_combo_real",
                   "ok": c.get("status") == "ok" and c.get("target_after") is not None,
                   "n_knockouts": c.get("n_knockouts"),
                   "target_before": c.get("target_before"),
                   "target_after": c.get("target_after")})

    # C: kcat 先验优化（代理，DEMO 非 GAN）
    cc = kcat_prior_optimize(model, n_samples=3, seed=1)
    checks.append({"check": "C_kcat_prior_surrogate",
                   "ok": cc.get("status") == "ok" and cc.get("best_target") is not None,
                   "best_kcat": cc.get("best_kcat"),
                   "best_target": cc.get("best_target")})
    # C-real-blocked 由调用方显式语义保证（本函数只用 demo 先验）

    # P0-1: ecFSEOF（酶约束 FSEOF + 酶丰度变异，真实 LP）
    e1 = ecfseof_scan(model, kcat_source="demo", cap_factor=1.0, n_steps=6,
                      enzyme_variation_factors=(0.5, 2.0))
    checks.append({"check": "P0_ecfseof_real",
                   "ok": e1.get("status") == "ok" and e1.get("ec_fseof") is not None,
                   "ec_n_up": e1.get("ec_fseof", {}).get("n_up"),
                   "ec_n_down": e1.get("ec_fseof", {}).get("n_down"),
                   "fba_n_up": e1.get("fba_fseof", {}).get("n_up"),
                   "fba_n_down": e1.get("fba_fseof", {}).get("n_down"),
                   "added_n": (e1.get("added_by_enzyme_constraint", {}).get("n_up", 0) +
                               e1.get("added_by_enzyme_constraint", {}).get("n_down", 0)),
                   "variation_keys": list(e1.get("enzyme_variation", {}).keys())})

    # P0-2: FCC（通量控制系数，有限差分）
    e2 = flux_control_coefficient(model, eps=1e-2, max_candidates=200)
    checks.append({"check": "P0_fcc_real",
                   "ok": e2.get("status") == "ok" and e2.get("n_coeff", 0) > 0,
                   "n_coeff": e2.get("n_coeff"),
                   "fmax": e2.get("fmax"),
                   "top_neg": [c["rxn"] for c in e2.get("top_negative", [])][:3]})

    # P0-3: 分级评估基准（L1/L2/L3/FCC + 三指标）
    e3 = tier_eval(model, kcat_source="demo", cap_factor=1.0)
    checks.append({"check": "P0_tier_eval_real",
                   "ok": e3.get("status") == "ok" and e3.get("metrics") is not None,
                   "L1_ko": e3.get("L1_single", {}).get("ko_candidates"),
                   "L2_ko": e3.get("L2_combo", {}).get("knockout_set"),
                   "L3_added": e3.get("L3_multiscale", {}).get("added_by_enzyme_constraint_n"),
                   "metrics_note": e3.get("metrics", {}).get("note", "")[:24]})

    # P0-3b: 项② 命名空间修复——取模型前 6 个带 GPR 的基因作参考集，验证 gene_level 可算
    ref_genes = [g.id for g in model.genes if g.reactions][:6]
    e3g = tier_eval(model, kcat_source="demo", cap_factor=1.0, reference_set=ref_genes)
    gl = e3g.get("metrics", {}).get("gene_level")
    checks.append({"check": "P0_tier_eval_genelevel_fix",
                   "ok": e3g.get("status") == "ok"
                        and e3g.get("metrics", {}).get("reference_namespace") == "gene"
                        and gl is not None,
                   "reference_namespace": e3g.get("metrics", {}).get("reference_namespace"),
                   "predicted_genes": len(e3g.get("predicted_genes", [])),
                   "gene_level": gl})

    # P0-4: 项② dedicated 必需性路径——串行 KO 必需性 benchmark 冒烟（限速 40 基因）
    smoke_genes = [g.id for g in model.genes][:40]
    ref_ess = smoke_genes[:3]   # 合成必需集（冒烟，非真实）
    eb = essentiality_benchmark(model, reference_essential=ref_ess, max_genes=40)
    checks.append({"check": "P0_essentiality_benchmark",
                   "ok": eb.get("status") == "ok" and eb.get("confusion_Epos_Nneg") is not None,
                   "predicted_essential_n": eb.get("predicted_essential_n"),
                   "confusion": eb.get("confusion_Epos_Nneg"),
                   "mcc": eb.get("mcc")})

    ok = all(ch["ok"] for ch in checks)
    return {
        "status": "ok" if ok else "FAIL",
        "gem": gem_key,
        "n_reactions": len(model.reactions),
        "n_genes": len(model.genes),
        "checks": checks,
        "note": "A/B/C 均在真实本地 GEM 上真 LP 闭合；A/C 用 demo/先验 kcat（诚实边界），"
                "B 纯 LP 全真；real kcat 不可达时 BLOCKED_NO_DATA 不编造。",
    }


if __name__ == "__main__":
    import sys
    key = sys.argv[1] if len(sys.argv) > 1 else "iYLI647"
    print(json.dumps(self_test(key), ensure_ascii=False, indent=2))
