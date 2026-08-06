# Xu lab · Mouse VPL TI Simulator v4

本包提供两个明确分开的计算层级：

- **Allen CCF 解析预览**：真实 100 µm 脑掩膜与左侧 VPL 标签、17,999 面透明脑表面、清晰三正交切面、四电极和轨迹、随电流与位置实时重算的 mTI 电场。用于坐标/参数筛查。
- **Registered EFMouse FEM**：在多组织四面体网格上分别求解 E1–E2 与 E3–E4 通道，再计算最大方向 TI 包络、定向包络、VPL/靶点球/全灰质统计与数值 QC。用于正式模板级计算。

默认立体定位约定是 AP 前方正、ML 右侧正、DV 背侧正，单位均为相对 bregma 的 mm。靶点：

```text
AP = -1.94, ML = -1.75, DV = -3.10 mm
```

## 1. Xu lab 默认四电极

| 电极 | 方向 | AP | ML | DV | 配对 |
|---|---:|---:|---:|---:|---|
| E1 | 前侧 | +0.14 | -1.75 | +0.50 | E1–E2 / AP |
| E2 | 后侧 | -4.02 | -1.75 | +0.50 | E1–E2 / AP |
| E3 | 内侧 | -1.94 | +0.33 | +0.50 | E3–E4 / ML |
| E4 | 外侧 | -1.94 | -3.83 | +0.50 | E3–E4 / ML |

每对电极的 AP/ML 中点均为 VPL 投影中心。默认半间距 2.08 mm；从 DV +0.50 mm 到靶点 DV -3.10 mm 的垂直落差是 3.60 mm，直线长度 `sqrt(2.08²+3.60²)=4.16 mm`，相对颅骨平面角度约 60.0°。E3 的 ML +0.33 mm 跨过中线；如果不允许跨中线，就必须旋转布局或改变角度/对称性。

网页把这四个点作为**请求的颅骨外表面坐标**显示。正式 FEM 会将各中心投影到已注册的真实外表面，并在 JSON 中报告每个投影残差；残差超过阈值时停止，而不是悄悄移动电极。

## 2. 修复 Python 3.12 / NumPy 安装失败

报错 `pkgutil.ImpImporter` 来自 Python 3.12 与旧 AllenSDK/旧 setuptools/旧 NumPy 构建链不兼容。v4 不依赖 AllenSDK；`nrrd` 模块由现代 `pynrrd` 提供。不要污染现有 `cxr_torch5090` 环境，建立独立 Python 3.11 环境：

```bash
conda deactivate
cd /mnt/e/ti/Xu_Lab_TI_VPL_Simulator_v4
bash create_environment.sh
conda activate ti_vpl_fem311
python -c "import numpy, scipy, nrrd, meshio, streamlit; print('environment OK')"
```

脚本会创建或更新 `ti_vpl_fem311`，不会删除其它环境。若 WSL 中 gmsh 缺少 OpenGL 动态库：

```bash
sudo apt-get update
sudo apt-get install -y libglu1-mesa libgl1
```

## 3. 启动本地 Web 软件

```bash
conda activate ti_vpl_fem311
cd /mnt/e/ti/Xu_Lab_TI_VPL_Simulator_v4
bash run_web.sh
```

浏览器打开 `http://localhost:8504`。只有点击“运行并提交新参数”后，草稿参数才成为新的 Run。可分别调整两对电极电流、`f1/f2`、十字半间距和旋转角；四电极坐标、轨迹、3D 场节点、切面与指标会重新计算。色标固定为 0–0.35 V/m，防止自动缩放掩盖真实变化。

### Streamlit Community Cloud 部署

仓库根目录已经包含云端入口和最小依赖：

```text
streamlit_app.py
requirements.txt
.streamlit/config.toml
```

在 [Streamlit Community Cloud](https://share.streamlit.io/) 中选择：

```text
Repository:  您上传后的 GitHub 仓库
Branch:      main
Main file:   streamlit_app.py
Python:      3.11（Advanced settings）
```

云端版运行 Allen 解剖约束解析预览。约 200 MB 的 EFMouse MAT/VTU 网格被 `.gitignore` 排除，避免超过 GitHub 文件限制；正式 registered FEM 应在本地/HPC 下载网格、完成配准后运行。不要在 Streamlit Cloud 上把未注册解析预览当成正式动物 FEM。

图中红/粉色实体是 **Allen 左侧 VPL 解剖标签**，不是电场焦点。彩色 3D 节点、切面热图和 0.045/0.10/0.20 V/m 等值线才是计算电场。默认通用 CCF 映射下，输入靶点距左侧 VPL 最近体素约 0.64 mm；软件保留并报告该不一致，不把靶点球冒充成 VPL。

## 4. 电流和频差怎样影响结果

- 在所用线性导电模型中，单通道场随对应注入电流线性变化；两通道不平衡会降低 `2·min(...)` 受限的一侧，并改变空间包络与质心。
- `Δf=|f2-f1|` 决定时间包络的拍频速度。对于当前频率不变的实数电导率、电准静态模型，只要 `Δf>0`，改变 Δf **不会改变空间包络形状或幅值**；`f1=f2` 时不存在 TI 拍频包络。
- 若要让载频改变空间场，必须使用经过测量/文献验证的复数阻抗与频散参数，分别求解各频率的复数场。不能仅凭频率数字人为改变热图。
- 对称几何与相等电流只保证几何中点和通道幅度平衡，不保证 VPL 是全脑最大值。应联合检查 VPL/灰质 P95 焦点性、高场质心到靶点距离、颅骨/表面热点与稳健性。

## 5. 下载并转换官方 EFMouse 网格

```bash
python scripts/download_efmouse_assets.py
python scripts/convert_efmouse_mesh.py
```

脚本下载官方 `Alekseichuk_mesh_clean.mat`、记录 SHA-256，并转换为 VTU。`pynrrd` 已在环境中，因此不会再出现 `ModuleNotFoundError: nrrd`。本压缩包不重复分发约 94 MB 源网格。

## 6. 建立 bregma → FEM 配准

```bash
cp registration/landmarks_example.csv registration/landmarks.csv
python scripts/fit_registration.py \
  --landmarks registration/landmarks.csv \
  --output-affine registration/bregma_to_fem_affine.txt \
  --output-qc registration/registration_qc.json \
  --maximum-rms-mm 0.20
```

CSV 列：`bregma_AP,bregma_ML,bregma_DV,fem_x,fem_y,fem_z`。至少 4 个非共面点，建议 6–10 个覆盖前后、左右和深度。示例恒等坐标仅演示格式，不能用于正式结论。正式计算需要个体影像或经验证模板的标志点、配准残差和切面 QC。

## 7. 运行高级 FEM

```bash
python scripts/run_advanced_fem.py --config config_advanced.yaml
```

正式模式会在缺少配准、配准 RMS 超限、电极投影距离超限、电流不守恒或注册 VPL 无灰质单元时停止。仅调试数据通路可显式使用：

```bash
python scripts/run_advanced_fem.py --config config_advanced.yaml --demo-unregistered
```

此输出会标记 `DEMO_UNREGISTERED_NOT_FOR_INFERENCE`。主要输出：

```text
outputs_advanced/ti_vpl_efmouse_registered.vtu
outputs_advanced/ti_vpl_efmouse_registered_metrics.json
```

模型求解 `∇·(σ∇φ)=0, E=-∇φ`，接触面采用严格配平的 Neumann 电流。VTU 包含两个载频矢量场、最大方向 TI 包络、AP/ML/DV 定向包络、载频 RMS、场夹角、Allen VPL 与靶点球。JSON 包含配准、电极投影、每通道电流、频率解释、电流守恒、线性残差、ROI 统计及焦点性。

## 8. 科学边界

这是模板级研究软件，不是动物剂量、疗效或安全性建议。发表或实验规划还需报告网格收敛、实测接触面积/阻抗、组织电导率敏感性、配准误差、求解残差、电流守恒、ROI 插值方式和独立软件/实测验证。建议引用 EFMouse 原论文（Sanchez-Romero et al., PLOS Computational Biology, 2025）与 Allen Mouse Brain CCF 的具体版本。
