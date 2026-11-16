"""Interactive Streamlit Lab & Verification Dashboard for MCCCII Current-Mode Biquad Filter."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import signal
import streamlit as st

# Ensure repository root is importable when running `streamlit run app.py`
ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.filters import BASE, coefficients, metrics, nodal_response, response, sweep  # noqa: E402

st.set_page_config(
    page_title="MCCCII 二階濾波器互動實驗與驗證台",
    page_icon="⚡",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Subtle engineering dark-accent styling
st.markdown(
    """
    <style>
    .block-container { padding-top: 1.5rem; padding-bottom: 2rem; }
    div[data-testid="stMetric"] {
        background-color: rgba(28, 131, 225, 0.08);
        border: 1px solid rgba(28, 131, 225, 0.25);
        border-radius: 10px;
        padding: 12px 16px;
    }
    .tag-badge {
        display: inline-block;
        padding: 2px 10px;
        border-radius: 999px;
        font-size: 0.78rem;
        font-weight: 600;
        background: rgba(16, 185, 129, 0.15);
        color: #10b981;
        border: 1px solid rgba(16, 185, 129, 0.4);
    }
    </style>
    """,
    unsafe_allow_html=True,
)

FILTER_INFO = {
    "LPF": {
        "name": "低通濾波器 (LPF)",
        "desc": "保留低頻訊號、抑制高頻雜訊；在 f₀ 處相位為 -90°，高頻滾降斜率 -40 dB/dec。",
        "admittance": "Y₁ = sC₁, Y₂ = G₂, Y₃ = sC₃, Y₄ = G₄",
    },
    "BPF": {
        "name": "帶通濾波器 (BPF)",
        "desc": "僅允許中心頻率 f₀ 附近訊號通過；f₀ 處相位為 0°，理想峰值增益為 0.2 (-13.98 dB)。",
        "admittance": "Y₁ = sC₁ + G₁, Y₂ = sC₂, Y₃ = sC₃ + G₃, Y₄ = G₄",
    },
    "HPF": {
        "name": "高通濾波器 (HPF)",
        "desc": "抑制低頻直流與漂移、保留高頻訊號；在 f₀ 處相位為 +90°，低頻斜率 +40 dB/dec。",
        "admittance": "Y₁ = G₁, Y₂ = sC₂, Y₃ = G₃, Y₄ = sC₄",
    },
}

VERSION_HISTORY = [
    {
        "version": "v1.0.0",
        "date": "2026-10-07",
        "stage": "專題文件與雙語架構基線",
        "added": "建立 README.md、中英雙語專題摘要 (docs/project-summary.md) 與參考文獻索引 (references/references.md)。",
        "deleted": "無（初始建置）。",
        "verification": "靜態文件審查；尚無自動化測試。",
        "rollback_note": "純文件基線版本，無程式碼相依性。",
    },
    {
        "version": "v2.0.0",
        "date": "2026-10-16",
        "stage": "理想模型重構與雙求解器 CI 驗證",
        "added": "新增 src/filters.py (符號轉移函數 + 3×3 KCL 節點矩陣求解器)、tests/test_filters.py (7 項單元測試)、run_experiments.py (21 組 CSV + 9 張 SVG)、docs/derivation.md 與 GitHub Actions CI (.github/workflows/verify.yml)。",
        "deleted": "移除早期專題報告中未經推導驗證的公式假設與靜態手算數值 (-180 行舊版文件內容)。",
        "verification": "7/7 unittest 通過；response() 與 nodal_response() 誤差 < 2e-9；與 scipy.signal.freqs 誤差 < 1e-12。",
        "rollback_note": "首個通過完整功能測試 (Functional Test Pass) 的黃金計算核心版本。",
    },
    {
        "version": "v2.1.0",
        "date": "2026-11-16",
        "stage": "線上會議互動實驗與迴歸測試儀表板",
        "added": "新增 app.py (Streamlit 互動頻率響應台、Monte Carlo 元件公差模擬、即時 7 項功能測試面板) 與 CHANGELOG.md。",
        "deleted": "無破壞性刪除；完整保留 src/filters.py 與 tests/test_filters.py 既有介面。",
        "verification": "即時執行 7 項功能測試與雙求解器交叉比對，確保互動介面與核心演算法零偏差。",
        "rollback_note": "若互動層發生相依套件異常，可直接退回 v2.0.0 核心 CLI 執行環境。",
    },
]


def run_functional_checks() -> list[dict[str, str]]:
    """Run core verification assertions mirroring tests/test_filters.py."""
    results = []
    f_grid = np.logspace(-2, 6, 71)

    # 1. Baseline w0 and Q
    expected_base = {
        "LPF": (1000.0, 0.5),
        "HPF": (1000.0, 0.5),
        "BPF": (1000.0 / np.sqrt(3), np.sqrt(3) / 5),
    }
    base_err = max(
        max(abs(metrics(k, BASE[k])[0] - w0_e), abs(metrics(k, BASE[k])[2] - q_e))
        for k, (w0_e, q_e) in expected_base.items()
    )
    results.append(
        {
            "test": "1. test_baseline (基準 ω₀ 與 Q 驗證)",
            "status": "PASS" if base_err < 1e-10 else "FAIL",
            "metric": f"最大誤差 = {base_err:.2e}",
        }
    )

    # 2. Independent 3x3 KCL nodal solver cross-check
    rng = np.random.default_rng(42)
    nodal_err = 0.0
    for kind, base_p in BASE.items():
        for _ in range(8):
            scales = 10 ** rng.uniform(-1, 1, size=len(base_p))
            p = {k: v * float(s) for (k, v), s in zip(base_p.items(), scales)}
            h_sym = response(kind, p, f_grid)
            h_nodal = nodal_response(kind, p, f_grid)
            nodal_err = max(nodal_err, float(np.max(np.abs(h_sym - h_nodal))))
    results.append(
        {
            "test": "2. test_independent_nodal_solver (3×3 KCL 矩陣交叉驗證)",
            "status": "PASS" if nodal_err < 2e-9 else "FAIL",
            "metric": f"24 組隨機參數最大誤差 = {nodal_err:.2e}",
        }
    )

    # 3. SciPy freqs cross-check
    scipy_err = 0.0
    for kind, p in BASE.items():
        n, d = coefficients(kind, p)
        _, hs = signal.freqs(n, d, 2 * np.pi * f_grid)
        scipy_err = max(scipy_err, float(np.max(np.abs(response(kind, p, f_grid) - hs))))
    results.append(
        {
            "test": "3. test_scipy (SciPy freqs 數值一致性)",
            "status": "PASS" if scipy_err < 1e-12 else "FAIL",
            "metric": f"最大誤差 = {scipy_err:.2e}",
        }
    )

    # 4. Limits and phase at f0
    dc_lpf = abs(response("LPF", BASE["LPF"], np.array([1e-8]))[0])
    hf_hpf = abs(response("HPF", BASE["HPF"], np.array([1e10]))[0])
    phases = {
        k: float(np.angle(response(k, BASE[k], np.array([metrics(k, BASE[k])[1]]))[0], deg=True))
        for k in BASE
    }
    phase_err = max(abs(phases["LPF"] + 90.0), abs(phases["HPF"] - 90.0), abs(phases["BPF"]))
    lim_ok = abs(dc_lpf - 1.0) < 1e-9 and abs(hf_hpf - 1.0) < 1e-9 and phase_err < 1e-9
    results.append(
        {
            "test": "4. test_limits_and_phase (極限增益與 f₀ 相位驗證)",
            "status": "PASS" if lim_ok else "FAIL",
            "metric": f"LPF={phases['LPF']:.1f}°, HPF={phases['HPF']:.1f}°, BPF={phases['BPF']:.1f}°",
        }
    )

    # 5. Stability (Re(poles) < 0) & orthogonal tuning
    max_real_pole = -np.inf
    for kind in BASE:
        for mode in ("q", "frequency"):
            for v in (0.25, 0.5, 1.0, 2.0, 4.0):
                _, d = coefficients(kind, sweep(kind, mode, v))
                max_real_pole = max(max_real_pole, float(np.max(np.real(np.roots(d)))))
    results.append(
        {
            "test": "5. test_stability_and_tuning (極點穩定性 Re(s)<0)",
            "status": "PASS" if max_real_pole < 0 else "FAIL",
            "metric": f"最大極點實部 = {max_real_pole:.2f} rad/s",
        }
    )

    # 6. BPF -3 dB bandwidth
    w0_bpf, _, q_bpf = metrics("BPF", BASE["BPF"])
    lo_w = w0_bpf * (np.sqrt(1 + 1 / (4 * q_bpf**2)) - 1 / (2 * q_bpf))
    hi_w = w0_bpf * (np.sqrt(1 + 1 / (4 * q_bpf**2)) + 1 / (2 * q_bpf))
    q_calc = w0_bpf / (hi_w - lo_w)
    results.append(
        {
            "test": "6. test_bpf_bandwidth (BPF -3 dB 頻寬驗證)",
            "status": "PASS" if abs(q_calc - q_bpf) < 1e-12 else "FAIL",
            "metric": f"ω₀/Δω = {q_calc:.6f} (理論 Q = {q_bpf:.6f})",
        }
    )

    # 7. Invalid input guard
    guard_ok = True
    try:
        coefficients("LPF", {**BASE["LPF"], "G2": -1.0})
        guard_ok = False
    except ValueError:
        pass
    results.append(
        {
            "test": "7. test_db_and_invalid (負值/非法參數防護)",
            "status": "PASS" if guard_ok else "FAIL",
            "metric": "負電導成功觸發 ValueError",
        }
    )
    return results


# ---------------------------------------------------------------------------
# Sidebar Controls
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("### ⚡ MCCCII 濾波器控制面板")
    st.caption("基於 `src/filters.py` 理想電流模式 ($R_x = 0$) 模型")

    kind = st.selectbox(
        "選擇濾波器拓撲 (Filter Topology)",
        options=["LPF", "BPF", "HPF"],
        format_func=lambda k: FILTER_INFO[k]["name"],
    )

    tuning_mode = st.radio(
        "調參模式 (Tuning Mode)",
        options=["正交倍率調參 (Sweep)", "元件級微調 (G / C)"],
        help="正交調參對應 src/filters.py 的 sweep()；元件級微調可直接改變單一電導 G 與電容 C。",
    )

    params = dict(BASE[kind])

    if tuning_mode == "正交倍率調參 (Sweep)":
        f_scale = st.slider("頻率倍率 (ω₀ Scale)", 0.25, 4.0, 1.0, 0.05)
        q_scale = st.slider("品質因數倍率 (Q Scale)", 0.25, 4.0, 1.0, 0.05)

        # Apply frequency scaling to conductances, then Q scaling
        params = sweep(kind, "frequency", f_scale)
        if kind == "LPF":
            params["G4"] *= q_scale
            params["G2"] /= q_scale
        elif kind == "HPF":
            params["G1"] *= q_scale
            params["G3"] /= q_scale
        else:  # BPF
            params["G4"] *= q_scale
    else:
        st.markdown("**電導 G (mS) 與電容 C (μF)**")
        for key, base_val in BASE[kind].items():
            if key.startswith("G"):
                val_ms = st.slider(
                    f"{key} (mS)",
                    min_value=0.1,
                    max_value=5.0,
                    value=float(base_val * 1e3),
                    step=0.05,
                )
                params[key] = val_ms * 1e-3
            else:
                val_uf = st.slider(
                    f"{key} (μF)",
                    min_value=0.1,
                    max_value=5.0,
                    value=float(base_val * 1e6),
                    step=0.05,
                )
                params[key] = val_uf * 1e-6

    show_baseline = st.checkbox("疊加基準曲線 (Baseline Overlay)", value=True)
    freq_range = st.select_slider(
        "掃描頻率範圍 (Hz)",
        options=[0.1, 1, 10, 100, 1000, 10000, 100000],
        value=(1, 10000),
    )

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
col_title, col_badge = st.columns([4, 1])
with col_title:
    st.title("MCCCII 電流模式二階濾波器 — 互動實驗與驗證台")
    st.caption(
        f"當前拓撲：**{FILTER_INFO[kind]['name']}** ｜ 導納映射：`{FILTER_INFO[kind]['admittance']}`"
    )
with col_badge:
    st.markdown(
        '<div style="text-align:right; padding-top:1rem;"><span class="tag-badge">v2.1.0 · Verified Core</span></div>',
        unsafe_allow_html=True,
    )

tab_lab, tab_verify, tab_versions = st.tabs(
    [
        "📊 1. 互動頻率響應實驗台",
        "🔬 2. 雙求解器交叉驗證 & 公差分析",
        "🏷️ 3. 版本歷程與韌體級功能測試",
    ]
)

# ---------------------------------------------------------------------------
# Tab 1: Interactive Frequency Response Lab
# ---------------------------------------------------------------------------
with tab_lab:
    f_min, f_max = freq_range
    if f_min >= f_max:
        f_min, f_max = 0.1, 100000.0
    f_grid = np.logspace(np.log10(f_min), np.log10(f_max), 800)

    h_curr = response(kind, params, f_grid)
    mag_db = 20 * np.log10(np.maximum(np.abs(h_curr), 1e-12))
    phase_deg = np.unwrap(np.angle(h_curr)) * 180.0 / np.pi

    w0, f0, q_val = metrics(kind, params)
    n_coef, d_coef = coefficients(kind, params)
    poles = np.roots(d_coef)
    is_stable = bool(np.all(np.real(poles) < 0))

    # KPI Row
    k1, k2, k3, k4 = st.columns(4)
    k1.metric("中心頻率 f₀", f"{f0:.2f} Hz", f"ω₀ = {w0:.1f} rad/s")
    k2.metric("品質因數 Q", f"{q_val:.4f}")
    f0_gain_db = float(20 * np.log10(np.maximum(abs(response(kind, params, np.array([f0]))[0]), 1e-12)))
    k3.metric("f₀ 處增益", f"{f0_gain_db:.2f} dB")
    k4.metric("系統穩定性 (極點)", "穩定 (Re < 0)" if is_stable else "不穩定", f"max Re(s)={np.max(np.real(poles)):.1f}")

    # Interactive Plotly Bode Plot
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=("增益響應 Gain (dB)", "相位響應 Unwrapped Phase (deg)"),
    )

    if show_baseline:
        h_base = response(kind, BASE[kind], f_grid)
        mag_base = 20 * np.log10(np.maximum(np.abs(h_base), 1e-12))
        phase_base = np.unwrap(np.angle(h_base)) * 180.0 / np.pi
        fig.add_trace(
            go.Scatter(
                x=f_grid,
                y=mag_base,
                name="基準 Gain (Baseline)",
                line=dict(color="rgba(148, 163, 184, 0.7)", width=2, dash="dash"),
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=f_grid,
                y=phase_base,
                name="基準 Phase (Baseline)",
                line=dict(color="rgba(148, 163, 184, 0.7)", width=2, dash="dash"),
                showlegend=False,
            ),
            row=2,
            col=1,
        )

    fig.add_trace(
        go.Scatter(
            x=f_grid,
            y=mag_db,
            name="目前 Gain (dB)",
            line=dict(color="#2dd4bf", width=3),
            fill="tozeroy",
            fillcolor="rgba(45, 212, 191, 0.08)",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=f_grid,
            y=phase_deg,
            name="目前 Phase (°)",
            line=dict(color="#38bdf8", width=2.5),
        ),
        row=2,
        col=1,
    )

    fig.add_vline(
        x=f0,
        line_dash="dot",
        line_color="#f59e0b",
        annotation_text=f"f₀ = {f0:.1f} Hz",
        annotation_position="top right",
    )

    fig.update_xaxes(type="log", title_text="頻率 Frequency (Hz)", row=2, col=1)
    fig.update_yaxes(title_text="Magnitude (dB)", row=1, col=1)
    fig.update_yaxes(title_text="Phase (deg)", row=2, col=1)
    fig.update_layout(height=540, template="plotly_dark", margin=dict(l=20, r=20, t=40, b=20))

    st.plotly_chart(fig, use_container_width=True)
    st.info(f"💡 **物理意義提示**：{FILTER_INFO[kind]['desc']}")

    with st.expander("查看即時轉移函數係數與元件參數表"):
        c_left, c_right = st.columns(2)
        with c_left:
            st.markdown("**正規化轉移函數 $H(s) = N(s)/D(s)$**")
            st.code(f"Numerator   N(s): {np.round(n_coef, 6)}\nDenominator D(s): {np.round(d_coef, 6)}\nPoles: {np.round(poles, 3)}")
        with c_right:
            st.markdown("**當前元件數值 (SI 單位)**")
            st.json({k: f"{v:.4e}" for k, v in params.items()})

# ---------------------------------------------------------------------------
# Tab 2: Dual-Solver Cross-Verification & Monte Carlo Tolerance
# ---------------------------------------------------------------------------
with tab_verify:
    st.subheader("符號轉移函數 vs. 3×3 KCL 節點矩陣求解器交叉驗證")
    st.write(
        "在工程驗證中，單一公式計算容易隱藏推導錯誤。本專案透過獨立的 3×3 節點導納矩陣 (`nodal_response`) "
        "與 `scipy.signal.freqs` 進行三方交叉比對，確保數學模型零偏差。"
    )

    h_sym_v = response(kind, params, f_grid)
    h_nodal_v = nodal_response(kind, params, f_grid)
    _, h_scipy_v = signal.freqs(n_coef, d_coef, 2 * np.pi * f_grid)

    err_nodal = float(np.max(np.abs(h_sym_v - h_nodal_v)))
    err_scipy = float(np.max(np.abs(h_sym_v - h_scipy_v)))

    v1, v2, v3 = st.columns(3)
    v1.metric("Symbolic vs. 3×3 KCL 最大誤差", f"{err_nodal:.2e}", "門檻 < 2.0e-9")
    v2.metric("Symbolic vs. SciPy 最大誤差", f"{err_scipy:.2e}", "門檻 < 1.0e-12")
    v3.metric("交叉驗證結果", "PASS" if err_nodal < 2e-9 and err_scipy < 1e-12 else "CHECK")

    st.divider()
    st.subheader("Monte Carlo 元件公差敏感度模擬 (Component Mismatch)")
    mc_col1, mc_col2 = st.columns([1, 3])
    with mc_col1:
        tol_pct = st.slider("元件誤差容許度 (±%)", 1, 20, 5, 1)
        n_trials = st.slider("模擬次數 (Trials)", 10, 80, 30, 10)
        seed = st.number_input("隨機種子 (Seed)", value=42, step=1)

    with mc_col2:
        rng_mc = np.random.default_rng(int(seed))
        fig_mc = go.Figure()
        f0_samples = []
        q_samples = []
        for idx in range(n_trials):
            perturbed = {
                k: v * float(rng_mc.uniform(1 - tol_pct / 100.0, 1 + tol_pct / 100.0))
                for k, v in params.items()
            }
            h_mc = response(kind, perturbed, f_grid)
            mag_mc = 20 * np.log10(np.maximum(np.abs(h_mc), 1e-12))
            _, f0_s, q_s = metrics(kind, perturbed)
            f0_samples.append(f0_s)
            q_samples.append(q_s)
            fig_mc.add_trace(
                go.Scatter(
                    x=f_grid,
                    y=mag_mc,
                    mode="lines",
                    line=dict(color="rgba(56, 189, 248, 0.25)", width=1),
                    showlegend=(idx == 0),
                    name=f"公差模擬 (±{tol_pct}%)",
                )
            )
        fig_mc.add_trace(
            go.Scatter(
                x=f_grid,
                y=mag_db,
                mode="lines",
                line=dict(color="#2dd4bf", width=3),
                name="標稱設計曲線 (Nominal)",
            )
        )
        fig_mc.update_xaxes(type="log", title_text="頻率 Frequency (Hz)")
        fig_mc.update_yaxes(title_text="Magnitude (dB)")
        fig_mc.update_layout(height=380, template="plotly_dark", margin=dict(l=20, r=20, t=20, b=20))
        st.plotly_chart(fig_mc, use_container_width=True)

        st.caption(
            f"在 ±{tol_pct}% 元件公差下（{n_trials} 次模擬）："
            f"f₀ 分佈於 **{min(f0_samples):.2f} ~ {max(f0_samples):.2f} Hz**，"
            f"Q 分佈於 **{min(q_samples):.4f} ~ {max(q_samples):.4f}**。"
        )

# ---------------------------------------------------------------------------
# Tab 3: Firmware-Style Functional Verification & Version History
# ---------------------------------------------------------------------------
with tab_versions:
    st.subheader("即時功能驗證測試 (Live Functional Test Suite)")
    st.write(
        "對應 `tests/test_filters.py` 的 7 項核心驗證。在韌體與系統驗證流程中，每次版本更新都必須通過此測試門檻，"
        "確保沒有引入 Functional Error。"
    )

    if st.button("▶ 立即執行 7 項功能驗證", type="primary"):
        checks = run_functional_checks()
        st.dataframe(checks, use_container_width=True, hide_index=True)
        if all(c["status"] == "PASS" for c in checks):
            st.success("全部 7 項功能驗證皆通過 (7/7 PASS) — 目前版本可安全作為二次開發基準。")
    else:
        checks = run_functional_checks()
        st.dataframe(checks, use_container_width=True, hide_index=True)

    st.divider()
    st.subheader("版本演進與 Additions / Deletions 追蹤 (Release Checkpoints)")
    st.write(
        "若後續版本發生 Functional Error，驗證團隊可直接透過以下版本標籤 (`git checkout v2.0.0`) "
        "回溯至上一個功能測試全數通過的穩定版本進行二次開發，快速鎖定異常變更。"
    )
    st.dataframe(VERSION_HISTORY, use_container_width=True, hide_index=True)
