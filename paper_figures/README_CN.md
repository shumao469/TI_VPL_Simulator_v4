# Xu lab · VPL–TI 论文图件

## 本次实际完成

根据附图“透明脑表面、内部场分布、靶区轮廓”的视觉方向，重新生成六组科学图。使用您现有 v4 解析公式和真实 Allen CCF 脑/VPL 资产；没有复制附图像素，也没有人为绘制靶区热点。

完整下载包包含：六张 600 dpi PNG、六份矢量 PDF、六份文字可编辑 SVG、三个透明背景的三维图层、英文图注、计算方法、源数据、八项数值测试及复现代码。图宽为 182.88 mm（7.2英寸），主图文字约7–11 pt。三维网格为投影后的矢量多边形；SVG文字保留为文本。

## 六组图的用途

1. **Figure 01**：透明脑及侧面、冠状、背侧三个剖切视角，金色为左 VPL，十字为请求靶点。
2. **Figure 02**：分别在请求靶点和 VPL 质心处截取三正交切面；没有 VPL 的切面会如实标明。
3. **Figure 03**：四电极理想布局、两个独立通道场与 mTI。DV +0.50 mm 是理想平面，不是重建后的真实颅骨面。
4. **Figure 04**：总电流固定1 mA，比较 0.25:0.75、0.50:0.50、0.75:0.25 mA；下方报告VPL、请求点、全脑P95变化。
5. **Figure 05**：差频10、40、100 Hz的载频信号与时间包络；不虚构频差引起的空间焦点移动。
6. **Figure 06**：请求点与解剖 VPL 的对应关系，以及脑和 VPL 的模型体素分布。

## 能否直接作为正式 FEM 结果发表？

目前这套图应表述为 **analytic model visualization / 解析模型示意与结果**。现有资产没有经过验证的个体/模板配准，也没有本次实际求得的多组织 FEM 电场。图中保留这一标识，英文图注已经说明模型边界。

默认模型得到的最大值不在 VPL，而在较浅位置。VPL平均值约0.1724、请求点约0.3316、全脑最大约3.0794（均为继承经验缩放系数的预览 V/m）。因此主图统一画相对值，不能把这些未标定值写成实际动物剂量。

请求靶点与最近左 VPL 体素中心相距0.6357 mm；这是通用坐标变换下的几何距离，不是已测量的配准误差，也不是到连续表面的精确距离。未做自动吸附或靶区平移。

正式结果需要您提供或在本地计算经配准的 FEM 输出，再用相同视觉规范绘制。相关输出通常为 `ti_vpl_efmouse_registered.vtu` 和配套 metrics JSON，但仍须检查配准、数组轴顺序、左侧ROI、电流守恒及网格收敛。

## 在您的 WSL 复现

```bash
cd /mnt/e/ti/TI_VPL_Simulator_v4
python3 -m venv .venv-figures
source .venv-figures/bin/activate
python -m pip install -r paper_figures/requirements.txt
python paper_figures/source/render_figures.py
python -m unittest discover -s paper_figures/source -p 'test_*.py' -v
```

若使用独立下载包，进入解压后同时包含 `assets/`、`ti_vpl/` 和 `paper_figures/` 的目录即可，路径按实际位置修改。无需 GPU、OpenGL 或 AllenSDK。

GitHub 提供代码、说明、数值摘要和轻量预览；完整下载包提供全部高清/矢量图与数值数组。GitHub代码运行后也可重新生成完整文件。

图中文字应在 SVG/PDF 中编辑；如改电流、电极或模型，应重新运行代码，不能只移动颜色热点或修改色标数值。
