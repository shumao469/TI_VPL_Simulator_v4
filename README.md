# Mouse TI Simulator · 小鼠 VPL 时间干涉电刺激仿真

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tivplsimulatorv4-dp7r8p99mgktdkfyjkmhyf.streamlit.app/)
[![CI](https://github.com/shumao469/TI_VPL_Simulator_v4/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/shumao469/TI_VPL_Simulator_v4/actions/workflows/ci.yml)

Xu Lab 的研究原型：结合 **Allen Mouse CCFv3 2017 解剖标签**、两通道四电极布局和时间干涉（Temporal Interference，TI）包络计算，探索小鼠丘脑腹后外侧核（VPL）附近的电场分布。网页提供快速解析预览，本地另提供基于 EFMouse 网格的有限元（FEM）流程。

**直接使用：[打开在线网页](https://tivplsimulatorv4-dp7r8p99mgktdkfyjkmhyf.streamlit.app/)。无需安装 Python。** 若应用休眠，按页面提示唤醒并等待加载。在线服务可用性取决于 Streamlit 托管状态。

<img width="1424" alt="Mouse TI Simulator 界面示例" src="https://github.com/user-attachments/assets/09c81e7b-06b4-4ce6-a3a8-5f96415368ac" />

## 1. 基于 Allen 脑所的哪个数据？

实际网页资源的来源记录位于 [assets/allen_ccf_100um_web.json](assets/allen_ccf_100um_web.json)，生成流程见 [scripts/build_atlas_assets.py](scripts/build_atlas_assets.py)。

| 项目 | 当前实现 |
|---|---|
| 图谱 | Allen Mouse Brain Common Coordinate Framework，第三版 **CCFv3** |
| 解剖标注版本 | **2017 年标注（ccf_2017）** |
| 网页实际分辨率 | **100 µm，各向同性** |
| 标注文件 | [annotation_100.nrrd](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_100.nrrd) |
| 配套平均模板 | [average_template_100.nrrd](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_100.nrrd) |
| 靶区 | VPL，Allen structure ID **718**；网页提取左半球标签 |
| 用途 | 全脑显示掩膜、透明脑表面、左 VPL 表面和切面轮廓 |

CCFv3 是三维参考空间，平均模板由 1,675 个标本脑构建，见 [Allen 官方文档](https://brain-map.org/support/documentation)及 [Wang et al., 2020](https://doi.org/10.1016/j.cell.2020.04.007)。图谱提供解剖参照；电导率与电场不是 Allen 图谱直接测得。

配置文件中仍有 `atlas_name: allen_mouse_25um` 字段，但当前网页及高级脚本实际读取 **100 µm** 标注，不能据该字段声称当前使用 25 µm 数据。

### 坐标约定

界面使用相对 bregma 的毫米：**AP 向前正、ML 向右正、DV 向背侧正**。默认请求靶点为 `AP=-1.94, ML=-1.75, DV=-3.10 mm`。

网页通用 QC 映射为 `AP=5.4−CCF_AP`、`ML=CCF_ML−5.7`、`DV=−CCF_DV`，单位 mm，CCF 数组轴序 AP/DV/ML。此近似不能代替个体或经验证模板配准。资源记录的默认靶点至左 VPL 最近体素距离约 **0.64 mm**；侧栏此值来自预生成资源，**不会随当前靶点输入重新计算**。

## 2. 基于什么方法仿真？

| 层级 | 方法 | 用途 |
|---|---|---|
| 在线/本地网页预览 | 均匀导体中点电流源—汇解析场，结合 Allen 解剖显示，再计算 TI 包络 | 电极、坐标、电流平衡与参数的定性探索 |
| 本地高级 FEM | EFMouse 多组织四面体网格，两通道电准静态有限元求解 | 给定材料、网格与配准假设下的模板电场和数值 QC |

### 解析预览

实现见 [ti_vpl/analytic.py](ti_vpl/analytic.py) 和 [app_advanced.py](app_advanced.py)。每个源—汇电极对计算：

```text
E(r) = k I / (4πσ) · [(r−rs)/|r−rs|³ − (r−rt)/|r−rt|³]
```

网页取 `σ=0.275 S/m`、固定缩放因子 `k=0.025`，距离分母下限 **0.45 mm** 用于限制近电极奇异值。底层 PreviewParameters 默认因子是 0.028，网页显式覆盖为 0.025。仓库未提供该因子的实验标定依据，预览绝对幅值不能作为动物剂量依据。

**Allen 掩膜约束显示范围，不参与多组织导电边界求解。** 三维脑解剖来自图谱，不意味着网页热图已经进行真实头模型 FEM。

### 有限元计算

入口为 [scripts/run_advanced_fem.py](scripts/run_advanced_fem.py)，Python 求解器见 [ti_vpl/fem.py](ti_vpl/fem.py)。使用 [EFMouse 官方仓库](https://github.com/klabhub/EFMouse)的 Alekseichuk_mesh_clean.mat 网格，转换后进行线性四面体有限元计算；不是直接调用 EFMouse MATLAB 工具箱，也未据此声明结果等同。

```text
∇·(σ∇φ) = 0
E = −∇φ
```

分别求解 E1–E2 和 E3–E4。接触面施加等量注入/流出的 Neumann 电流，其余边界无外加通量，固定一个节点电位。小系统稀疏直接求解；大系统共轭梯度，优先 AMG 预条件器，缺少 AMG 时回退至对角预条件器。

当前组织电导率（S/m）：灰质 0.275、脑脊液 1.654、骨 0.010、皮肤 0.465、眼 0.500。这些是需针对实验核实的配置假设。高级入口采用实数、各向同性、频率不变电导率，不包含神经元响应、热效应或频散仿真。

### TI 包络和频率

[ti_vpl/envelope.py](ti_vpl/envelope.py) 对两个实矢量场计算：

```text
定向包络：m(n) = 2 min(|n·E1|, |n·E2|)，|n|=1
最大方向包络：mTI = max_n m(n)
拍频：Δf = |f2−f1|
```

最大方向包络在各位置自由选择方向，是上界指标，不等于实际神经纤维刺激强度或激活概率。默认 f1=2000 Hz、f2=2010 Hz，拍频 10 Hz。在当前模型下，非零频差改变时间包络速度，不改变空间幅值；两频相等时网页/解析预览将 TI 拍频包络置零；当前高级 FEM 入口未做同样的等频分支，正式 TI 分析应使用不同载频，不能把等频时输出的空间公式值解释为拍频。

## 3. 如何操作在线网页？

1. 打开在线链接，等待界面和脑模型加载。
2. 在左侧输入靶点 AP/ML/DV，调整两通道电流、f1/f2、十字电极半间距和旋转角。
3. 点击 **“运行并提交新参数”**，检查 Run 编号和电极坐标表；只修改输入不会提交新 Run。
4. 拖动三维图旋转、滚轮缩放，结合冠状、矢状、水平切面查看电场与解剖位置。
5. 查看 target mTI（靶点单点解析包络）、pair current balance（两通道较小电流/较大电流）和 beat frequency。电流平衡百分比不是数值电流守恒残差。

**红/粉色实体是左 VPL 解剖标签，不是电场焦点。** 彩色节点和切面热图才是计算场。色标固定 0–0.35 V/m，超限可能颜色饱和，不能由色标上限推断最大值。

默认每通道 0.50 mA，半间距 2.08 mm、旋转 0°、电极 DV=+0.50 mm：

| 电极 | AP (mm) | ML (mm) | DV (mm) | 通道 |
|---|---:|---:|---:|---|
| E1 | +0.14 | −1.75 | +0.50 | E1–E2 |
| E2 | −4.02 | −1.75 | +0.50 | E1–E2 |
| E3 | −1.94 | +0.33 | +0.50 | E3–E4 |
| E4 | −1.94 | −3.83 | +0.50 | E3–E4 |

E3 跨过中线；此布局仅为可调整初值。等电流和几何对称不保证 VPL 是全脑场最大值。

## 4. 如何下载并在本地打开？

### 只运行网页预览（Windows / Linux / macOS）

安装 Git 和 **Python 3.11**，在终端执行：

```bash
git clone https://github.com/shumao469/TI_VPL_Simulator_v4.git
cd TI_VPL_Simulator_v4
python -m venv .venv
```

激活环境：Windows PowerShell 用 `.\.venv\Scripts\Activate.ps1`；Linux/macOS/WSL 用 `source .venv/bin/activate`。然后：

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run streamlit_app.py --server.port 8504
```

浏览器打开 **http://localhost:8504**，终端保持运行；Ctrl+C 停止。已激活环境下也可运行 Windows 的 run_web.bat，或 Linux/WSL 的 `bash run_web.sh`。端口占用时改为 8505 并打开对应地址。预览使用仓库内解剖资源，无需 AllenSDK 或 EFMouse 大网格。

### 完整 FEM 环境（建议 Linux / WSL + Conda/Miniforge）

在仓库根目录执行：

```bash
bash create_environment.sh
conda activate ti_vpl_fem311
bash run_web.sh
```

脚本按 environment.yml 创建或更新 ti_vpl_fem311 并安装项目。若 gmsh/VTK 缺少 OpenGL 库，Ubuntu/WSL 可安装 libglu1-mesa 和 libgl1。完整环境较大，安装时间取决于网络与依赖解析。

## 5. 如何运行本地 FEM？

以下在仓库根目录和完整环境中执行。

**第一步：下载并转换网格。** 大型 MAT/VTU 未随仓库提供；脚本记录下载 SHA-256。

```bash
python scripts/download_efmouse_assets.py
python scripts/convert_efmouse_mesh.py
```

**第二步：准备真实标志点并配准。** registration/landmarks_example.csv 仅示例格式，不能直接充当实验配准。准备 registration/landmarks.csv，列为 bregma_AP,bregma_ML,bregma_DV,fem_x,fem_y,fem_z，两侧单位 mm。至少 4 个非共面对应点；建议更多空间分布充分的点，并独立验证。4 点仿射拟合低残差本身不证明泛化准确。

```bash
python scripts/fit_registration.py --landmarks registration/landmarks.csv --output-affine registration/bregma_to_fem_affine.txt --output-qc registration/registration_qc.json --maximum-rms-mm 0.20
```

**第三步：检查配置并求解。** 核对 config_advanced.yaml 的坐标、电流、频率、电导率、接触面和输出目录。

```bash
python scripts/run_advanced_fem.py --config config_advanced.yaml
```

网页下方检测到网格和配准文件后，也可点 **“运行完整注册 FEM”**，当前界面参数会写入 outputs_advanced/web_run_config.yaml 后求解。ready 仅表示文件存在，内容有效性需在计算阶段检查。

主要输出：

- outputs_advanced/ti_vpl_efmouse_registered.vtu：两通道电场、TI 包络和 ROI 单元数据，可用 ParaView 查看。
- outputs_advanced/ti_vpl_efmouse_registered_metrics.json：配准、投影、ROI、电流与线性残差等；网页可下载最近结果 JSON。

缺少配准、配准 RMS 超限、电极投影距离超限或 VPL 无灰质单元时，正式入口会停止。`--demo-unregistered` 仅调试通路，不用于正式推断。

**当前 ROI 限制：** 高级脚本仍经通用 bregma→CCF 映射采样 100 µm 标注，Allen_VPL 按 ID 718 取样，未像网页一样显式限制左半球。bregma→FEM 配准通过不代表个体 Allen 图谱配准完成，解读左 VPL 统计前需核验半球、轴序及 ROI 对齐。

## 6. 测试、故障排查与适用边界

```bash
python -m pip install -e . pytest
python -m pytest
```

测试覆盖默认坐标、电极几何、包络性质、电流线性缩放、单四面体求解与仿射恢复；CI 另做语法检查。**CI 通过表示安装与现有测试通过，不代表动物或真实头模型验证完成。**

| 问题 | 处理 |
|---|---|
| 参数修改后图未更新 | 点击“运行并提交新参数”，核对 Run 编号 |
| FEM 按钮不可用 | 下载/转换网格，并准备真实配准文件 |
| No module named nrrd | 在当前环境安装 pynrrd |
| Multiple top-level packages discovered | 更新仓库；pyproject.toml 已明确只发现 ti_vpl 包 |
| 云端无 FEM 网格 | 在线使用预览，本地配置完整流程 |

科研使用仍需网格收敛、接触/组织参数敏感性、配准/ROI 检查和独立软件或实测对照。本程序不提供神经激活、疗效或安全剂量结论。[中文补充说明](README_CN.md)包含更多操作背景，图谱和实现细节以本页及当前源码为准。

## 7. 数据与方法来源

- [Allen 官方 CCF 文档](https://brain-map.org/support/documentation)：CCFv3、2017 标注与模板。
- [Wang et al. (2020), Cell](https://doi.org/10.1016/j.cell.2020.04.007)：The Allen Mouse Brain Common Coordinate Framework: A 3D Reference Atlas。
- [EFMouse 官方仓库](https://github.com/klabhub/EFMouse)及 [Sanchez-Romero, Akyuz & Krekelberg (2025), PLOS Computational Biology](https://doi.org/10.1371/journal.pcbi.1013471)：网格和原工具箱背景。
- 本仓库解析预览、Python FEM 与包络分别见 ti_vpl/analytic.py、ti_vpl/fem.py、ti_vpl/envelope.py。使用外部资源请遵循其许可和引用要求。
