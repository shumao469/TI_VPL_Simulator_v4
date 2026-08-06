#!/usr/bin/env python3
"""Xu Lab Allen CCF preview and launcher for registered EFMouse FEM."""
from __future__ import annotations

import base64
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import plotly.graph_objects as go
import streamlit as st
import yaml

from ti_vpl.analytic import point_pair_field
from ti_vpl.envelope import maximum_orientation_envelope

ROOT = Path(__file__).resolve().parent
ASSET_PATH = ROOT / "assets/allen_ccf_100um_web.json"
CONFIG_PATH = ROOT / "config_advanced.yaml"
FIELD_CAP_V_PER_M = 0.35
FIELD_ISO_V_PER_M = (0.045, 0.10, 0.20)


@st.cache_resource
def load_atlas() -> dict:
    asset = json.loads(ASSET_PATH.read_text(encoding="utf-8"))
    count = int(np.prod(asset["shape_AP_DV_ML"]))
    for source, target in (("brain_mask_bitpack_little_base64", "brain"), ("left_vpl_mask_bitpack_little_base64", "vpl")):
        packed = np.frombuffer(base64.b64decode(asset[source]), dtype=np.uint8)
        asset[target] = np.unpackbits(packed, bitorder="little")[:count].reshape(asset["shape_AP_DV_ML"]).astype(bool)
    return asset


def electrode_centres(target: np.ndarray, spacing: float, rotation_deg: float, surface_dv: float = 0.50) -> np.ndarray:
    theta = np.deg2rad(rotation_deg)
    ap_axis = np.asarray([np.cos(theta), np.sin(theta)])
    ml_axis = np.asarray([-np.sin(theta), np.cos(theta)])
    center = target[:2]
    xy = np.vstack([center + spacing * ap_axis, center - spacing * ap_axis,
                    center + spacing * ml_axis, center - spacing * ml_axis])
    return np.column_stack([xy, np.full(4, surface_dv)])


def field(points: np.ndarray, target: np.ndarray, current_1_mA: float, current_2_mA: float,
          spacing: float, rotation_deg: float, f1: float, f2: float) -> np.ndarray:
    electrodes = electrode_centres(target, spacing, rotation_deg)
    common = {"conductivity_S_per_m": 0.275, "skull_transfer_factor": 0.025}
    e1 = point_pair_field(points, electrodes[0], electrodes[1], current_mA=current_1_mA, **common)
    e2 = point_pair_field(points, electrodes[2], electrodes[3], current_mA=current_2_mA, **common)
    modulation = maximum_orientation_envelope(e1, e2)
    return np.zeros_like(modulation) if np.isclose(f1, f2) else modulation


def sample_mask(asset: dict, points: np.ndarray, name: str) -> np.ndarray:
    p = np.asarray(points)
    ijk = np.rint(np.column_stack([(5.4 - p[:, 0]) / 0.1, -p[:, 2] / 0.1, (p[:, 1] + 5.7) / 0.1])).astype(int)
    shape = np.asarray(asset["shape_AP_DV_ML"])
    valid = np.all((ijk >= 0) & (ijk < shape), axis=1)
    result = np.zeros(len(p), dtype=bool)
    q = ijk[valid]
    result[valid] = asset[name][q[:, 0], q[:, 1], q[:, 2]]
    return result


def slice_figure(asset: dict, kind: str, target: np.ndarray, params: tuple[float, ...]) -> go.Figure:
    if kind == "冠状":
        x, y = np.linspace(-5.7, 5.7, 260), np.linspace(-8.0, 0.0, 220)
        xx, yy = np.meshgrid(x, y); points = np.column_stack([np.full(xx.size, target[0]), xx.ravel(), yy.ravel()])
        marker = (target[1], target[2]); labels = ("ML (mm)", "DV (mm)")
    elif kind == "矢状":
        x, y = np.linspace(-7.7, 5.4, 280), np.linspace(-8.0, 0.0, 220)
        xx, yy = np.meshgrid(x, y); points = np.column_stack([xx.ravel(), np.full(xx.size, target[1]), yy.ravel()])
        marker = (target[0], target[2]); labels = ("AP (mm)", "DV (mm)")
    else:
        x, y = np.linspace(-5.7, 5.7, 260), np.linspace(5.4, -7.7, 280)
        xx, yy = np.meshgrid(x, y); points = np.column_stack([yy.ravel(), xx.ravel(), np.full(xx.size, target[2])])
        marker = (target[1], target[0]); labels = ("ML (mm)", "AP (mm)")
    values = field(points, target, *params)
    brain = sample_mask(asset, points, "brain")
    vpl = sample_mask(asset, points, "vpl")
    z = np.where(brain, values, np.nan).reshape(xx.shape)
    vpl_z = np.where(vpl, 1.0, np.nan).reshape(xx.shape)
    fig = go.Figure(go.Heatmap(x=x, y=y, z=z, zmin=0, zmax=FIELD_CAP_V_PER_M,
                               zsmooth=False, colorscale="Turbo", colorbar={"title": "TI V/m", "thickness": 9}))
    for level in FIELD_ISO_V_PER_M:
        fig.add_trace(go.Contour(x=x, y=y, z=z, contours={"start": level, "end": level, "size": 1e-6, "coloring": "lines"},
                                 line={"color": "white", "width": 1}, showscale=False, hoverinfo="skip", name=f"{level:.3f} V/m"))
    fig.add_trace(go.Contour(x=x, y=y, z=vpl_z, contours={"start": .5, "end": .5, "size": .1, "coloring": "lines"},
                             line={"color": "#d33f58", "width": 2}, showscale=False, hoverinfo="skip", name="Allen VPL"))
    fig.add_trace(go.Scatter(x=[marker[0]], y=[marker[1]], mode="markers", marker={"size": 12, "color": "rgba(0,0,0,0)", "line": {"color": "white", "width": 3}}, name="target"))
    fig.update_layout(height=340, margin=dict(l=4, r=4, t=28, b=4), xaxis_title=labels[0], yaxis_title=labels[1],
                      title=kind, showlegend=False, paper_bgcolor="white", plot_bgcolor="#f4f7f8", font_color="#15263b")
    return fig


def mesh_trace(mesh: dict, color: str, opacity: float, name: str) -> go.Mesh3d:
    v, f = np.asarray(mesh["vertices"]), np.asarray(mesh["faces"], dtype=int)
    return go.Mesh3d(x=v[:, 1], y=v[:, 0], z=v[:, 2], i=f[:, 0], j=f[:, 1], k=f[:, 2],
                     color=color, opacity=opacity, flatshading=False, lighting={"ambient": .70, "diffuse": .65, "specular": .18},
                     lightposition={"x": 100, "y": 80, "z": 120}, name=name, hoverinfo="name")


def figure_3d(asset: dict, target: np.ndarray, params: tuple[float, ...]) -> go.Figure:
    ap, ml, dv = np.meshgrid(np.arange(-7.2, 4.91, .36), np.arange(-5.4, 5.41, .36), np.arange(-7.7, -.09, .36), indexing="ij")
    points = np.column_stack([ap.ravel(), ml.ravel(), dv.ravel()])
    inside = sample_mask(asset, points, "brain")
    points = points[inside]
    values = field(points, target, *params)
    hot = values >= FIELD_ISO_V_PER_M[0]
    electrodes = electrode_centres(target, params[2], params[3])
    fig = go.Figure([
        mesh_trace(asset["brain_mesh"], "#80c4cf", .16, "Allen whole brain (17,999 faces)"),
        mesh_trace(asset["left_vpl_mesh"], "#df4d61", .62, "Allen left VPL anatomy"),
        go.Scatter3d(x=points[hot, 1], y=points[hot, 0], z=points[hot, 2], mode="markers", name="calculated mTI ≥0.045 V/m",
                     marker={"size": 2.6, "color": values[hot], "cmin": 0, "cmax": FIELD_CAP_V_PER_M, "colorscale": "Turbo", "opacity": .55,
                             "colorbar": {"title": "mTI V/m", "len": .50, "thickness": 10}}),
        go.Scatter3d(x=[target[1]], y=[target[0]], z=[target[2]], mode="markers", name="VPL target",
                     marker={"size": 8, "color": "white", "symbol": "cross", "line": {"color": "#10243e", "width": 3}}),
        go.Scatter3d(x=electrodes[:, 1], y=electrodes[:, 0], z=electrodes[:, 2], mode="markers+text", text=["E1", "E2", "E3", "E4"],
                     textposition="top center", name="requested skull-surface electrodes",
                     marker={"size": 9, "color": ["#2474cf", "#2474cf", "#e25555", "#e25555"], "line": {"color": "white", "width": 1}}),
    ])
    for index, electrode in enumerate(electrodes):
        fig.add_trace(go.Scatter3d(x=[electrode[1], target[1]], y=[electrode[0], target[0]], z=[electrode[2], target[2]],
                                   mode="lines", line={"color": "#6c7e88", "width": 3, "dash": "dash"},
                                   name=f"E{index + 1}→VPL trajectory", showlegend=False, hoverinfo="skip"))
    fig.update_layout(height=610, margin=dict(l=0, r=0, t=34, b=0), paper_bgcolor="white", font_color="#15263b",
                      scene={"xaxis_title": "ML", "yaxis_title": "AP", "zaxis_title": "DV", "aspectmode": "data", "bgcolor": "white",
                             "xaxis": {"backgroundcolor": "white", "gridcolor": "#d9e3e7"},
                             "yaxis": {"backgroundcolor": "white", "gridcolor": "#d9e3e7"},
                             "zaxis": {"backgroundcolor": "white", "gridcolor": "#d9e3e7"}},
                      title="Xu lab · Allen CCF 透明 3D · VPL anatomy · calculated mTI field")
    return fig


def beat_figure(target_field: float, f1: float, f2: float) -> go.Figure:
    beat = abs(f2 - f1)
    duration = min(0.50, max(0.08, 3.0 / beat)) if beat > 0 else 0.20
    time = np.linspace(0, duration, 800)
    envelope = target_field * np.abs(np.cos(np.pi * beat * time)) if beat > 0 else np.zeros_like(time)
    fig = go.Figure(go.Scatter(x=time, y=envelope, line={"color": "#0b9186", "width": 2}, fill="tozeroy", fillcolor="rgba(11,145,134,.12)"))
    fig.update_layout(height=190, margin=dict(l=28, r=12, t=28, b=28), title=f"靶点时间包络 · Δf={beat:g} Hz",
                      xaxis_title="time (s)", yaxis_title="mTI (V/m)", paper_bgcolor="white", plot_bgcolor="white")
    return fig


def main() -> None:
    st.set_page_config(page_title="Xu lab | Mouse VPL TI Simulator", page_icon="🧠", layout="wide")
    asset = load_atlas()
    st.title("Xu lab · Mouse VPL TI Simulator")
    st.caption("Allen CCF anatomy-aware analytic preview + registered EFMouse FEM · AP前+/ML右+/DV背侧+ · mm from bregma")
    with st.sidebar:
        st.subheader("靶点与四电极")
        with st.form("parameters"):
            ap = st.number_input("target AP", -8.0, 5.0, -1.94, .01)
            ml = st.number_input("target ML", -6.0, 6.0, -1.75, .01)
            dv = st.number_input("target DV", -8.0, .5, -3.10, .01)
            current_1 = st.slider("E1–E2 current (mA)", .05, 1.0, .50, .05)
            current_2 = st.slider("E3–E4 current (mA)", .05, 1.0, .50, .05)
            f1 = st.number_input("f1 (Hz)", 100.0, 10000.0, 2000.0, 1.0)
            f2 = st.number_input("f2 (Hz)", 100.0, 10000.0, 2010.0, 1.0)
            spacing = st.slider("cross half spacing (mm)", 1.0, 3.0, 2.08, .01)
            rotation = st.slider("montage rotation (deg)", -45.0, 45.0, 0.0, 1.0)
            submitted = st.form_submit_button("运行并提交新参数", type="primary", width="stretch")
        if submitted or "params" not in st.session_state:
            st.session_state.params = (ap, ml, dv, current_1, current_2, spacing, rotation, f1, f2)
            st.session_state.run_id = st.session_state.get("run_id", 0) + 1
        p = st.session_state.params
        st.success(f"Run #{st.session_state.run_id} · target {p[0]:.2f}/{p[1]:.2f}/{p[2]:.2f}")
        target = np.asarray(p[:3], dtype=float)
        electrodes = electrode_centres(target, p[5], p[6])
        st.dataframe({"electrode": ["E1", "E2", "E3", "E4"], "AP": electrodes[:, 0], "ML": electrodes[:, 1], "DV": electrodes[:, 2]}, hide_index=True)
        drop = .50 - target[2]
        trajectory = np.hypot(drop, p[5])
        angle = np.degrees(np.arctan2(drop, p[5]))
        st.caption(f"vertical drop {drop:.2f} mm · trajectory {trajectory:.2f} mm · angle {angle:.1f}° from skull plane")
        st.metric("generic CCF target→left VPL residual", f"{asset['generic_mapping_target_to_vpl_distance_mm']:.2f} mm")

    target = np.asarray(p[:3], dtype=float)
    field_params = (p[3], p[4], p[5], p[6], p[7], p[8])
    st.plotly_chart(figure_3d(asset, target, field_params), width="stretch")
    target_field = float(field(target.reshape(1, 3), target, *field_params)[0])
    c1, c2, c3 = st.columns(3)
    c1.metric("target mTI", f"{target_field:.3f} V/m")
    c2.metric("pair current balance", f"{min(p[3], p[4]) / max(p[3], p[4]) * 100:.0f}%")
    c3.metric("beat frequency", f"{abs(p[8] - p[7]):g} Hz")
    st.plotly_chart(beat_figure(target_field, p[7], p[8]), width="stretch")
    st.info("红/粉色实体是 Allen 左侧 VPL 解剖标签，不是计算焦点。彩色节点与切面热图才是解析 mTI 电场；固定色标 0–0.35 V/m。Δf 只改变时间包络速度；在当前电阻性电准静态模型中，非零 Δf 不改变空间包络。f1=f2 时没有 TI 拍频包络。")
    columns = st.columns(3)
    for col, kind in zip(columns, ("冠状", "矢状", "水平")):
        with col:
            st.plotly_chart(slice_figure(asset, kind, target, field_params), width="stretch")

    st.subheader("Registered EFMouse full FEM")
    mesh_ok = (ROOT / "data/efmouse/Alekseichuk_mesh_clean.vtu").exists()
    registration_ok = (ROOT / "registration/bregma_to_fem_affine.txt").exists() and (ROOT / "registration/registration_qc.json").exists()
    q1, q2, q3 = st.columns(3)
    q1.metric("EFMouse mesh", "ready" if mesh_ok else "missing")
    q2.metric("registration QC", "ready" if registration_ok else "missing")
    q3.metric("mode", "registered FEM" if mesh_ok and registration_ok else "analytic preview only")
    st.warning("Equal currents and geometric symmetry only balance the two pair amplitudes; they do not prove the VPL is the whole-brain maximum. Inspect VPL/gray P95 focality, high-field centroid and surface hot spots in the registered FEM output.")
    if st.button("运行完整注册 FEM", disabled=not (mesh_ok and registration_ok), type="primary"):
        run_cfg = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
        run_cfg["coordinates"]["requested_target_mm"] = {"AP": float(target[0]), "ML": float(target[1]), "DV": float(target[2])}
        run_cfg["stimulation"]["channel_1_current_mA"] = float(p[3])
        run_cfg["stimulation"]["channel_2_current_mA"] = float(p[4])
        run_cfg["stimulation"]["f1_Hz"], run_cfg["stimulation"]["f2_Hz"] = float(p[7]), float(p[8])
        updated = electrode_centres(target, p[5], p[6])
        run_cfg["stimulation"]["channel_1"]["source_bregma_mm"] = updated[0].tolist()
        run_cfg["stimulation"]["channel_1"]["return_bregma_mm"] = updated[1].tolist()
        run_cfg["stimulation"]["channel_2"]["source_bregma_mm"] = updated[2].tolist()
        run_cfg["stimulation"]["channel_2"]["return_bregma_mm"] = updated[3].tolist()
        output_cfg = ROOT / "outputs_advanced/web_run_config.yaml"
        output_cfg.parent.mkdir(parents=True, exist_ok=True)
        output_cfg.write_text(yaml.safe_dump(run_cfg, sort_keys=False), encoding="utf-8")
        with st.spinner("registered FEM running …"):
            result = subprocess.run([sys.executable, str(ROOT / "scripts/run_advanced_fem.py"), "--config", str(output_cfg)], cwd=ROOT, text=True, capture_output=True)
        st.error(result.stderr or result.stdout) if result.returncode else st.success(result.stdout)
    metrics_path = ROOT / "outputs_advanced/ti_vpl_efmouse_registered_metrics.json"
    if metrics_path.exists():
        st.download_button("下载最近 FEM 指标 JSON", metrics_path.read_bytes(), "ti_vpl_efmouse_registered_metrics.json")
        st.json(json.loads(metrics_path.read_text(encoding="utf-8")), expanded=False)


if __name__ == "__main__":
    main()
