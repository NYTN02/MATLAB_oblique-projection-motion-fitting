# 非理想斜抛运动发射角预测系统

## 简介

本 MATLAB 代码通过拟合实验数据，得到非理想斜抛运动的发射角（方位角 θ）预测模型。模型包含物理解析式、k 函数、θ 扰动及分段校正，可显著降低预测误差。

## 文件结构

- `theta_correction.m` : 主程序
- `step6_results.mat` : k 函数优化参数（运行后生成）
- `step7_results.mat` : θ 扰动参数（运行后生成）
- `step8_results.mat` : 分段校正参数（运行后生成）

## 环境要求

- MATLAB R2016b+
- Curve Fitting Toolbox
- Optimization Toolbox
- Global Optimization Toolbox（可选，用于 patternsearch）

## 快速开始

将 `theta_correction.m` 放入 MATLAB 工作目录，在命令窗口输入 `theta_correction` 运行。约 1-2 分钟后生成三个 `.mat` 参数文件。

## 算法步骤

- 1. **常数 k 拟合**  
基于解析式
```
θ_theory = atan2( 1 ± √(1 - k·r·(2·sinα + k·r·cos²α)) ,  k·r·cosα ) × 180/π
```
搜索最优 k。

- 2. **k 函数优化（步骤6）**  
将 k 扩展为 r, α 的函数：  
```
k = k_base * (1 + a·r + b·α + c·r·α + d·r² + e·α²)
```
输出参数 `k_base, a, b, c, d, e` 保存于 `step6_results.mat`。

- 3. **θ 扰动优化（步骤7）**  
对理论 θ 添加扰动：  
```
θ_pert = θ_theory + A1 + A2·r + A3·sinα + A4·sin2α + A5·r·sinα
```
输出参数 `A1~A5` 保存于 `step7_results.mat`。

- 4. **分段校正（步骤8）**  
根据 θ_pert 是否 ≤45° 分别采用两组线性模型（含 θ、θ²、r、r²、sinα、cosα 等项）计算校正量 Δθ。  
输出两套系数保存于 `step8_results.mat`。

- 5.**对比多项式拟合**

## 最终预测公式
将运行得到的参数代入以下表达式，即可计算任意目标(r, α) 对应的 θ：
- **1. 原始 θ 解析式（取小角度解）**  
```θ_theory = atan2( 1 - √(1 ± k·r·(2·sinα + k·r·cos²α)) ,  k·r·cosα ) × 180/π```
- **2. k 函数表达式**（参数来自步骤6，保存于 step6_results.mat）  
```k = k_base × (1 + a·r + b·α + c·r·α + d·r² + e·α²)```
- **3. θ 扰动表达式**（参数来自步骤7，保存于 step7_results.mat）  
```θ_pert = θ_theory + A1 + A2·r + A3·sinα + A4·sin(2α) + A5·r·sinα```
- **4. θ 校正函数表达式**（参数来自步骤8，保存于 step8_results.mat）  
- 若 θ_pert ≤ 45°：  
```
Δθ = c0 + c1·θ_pert + c2·θ_pert² + c3·(r/300) + c4·(r/300)² + c5·sinα + c6·cosα + c7·sin2α + c8·(r/300)·sinα + c9·(α/90)
```
- 若 θ_pert > 45°：  
```
Δθ = d0 + d1·θ_pert + d2·θ_pert² + d3·(r/300) + d4·(r/300)² + d5·sinα + d6·sin2α
```

**最终预测值**：  
```
θ_pred = θ_pert + Δθ
```## 输出图形说明

程序运行过程中会依次弹出 8 个图形窗口，各图说明如下：

| 图号 | 名称 | 说明 |
| :---: | :--- | :--- |
| 1 | 多项式拟合曲面 | 三维数据点与拟合曲面，用于观察整体趋势 |
| 2 | 大范围 k 搜索曲线 | RMSE 随 k 变化，初步定位最优 k 区间 |
| 3 | 解析式拟合曲面 | 两个理论解（红/蓝）与实验点（绿）对比 |
| 4 | 精细 k 搜索曲线 | 局部放大 RMSE 曲线，精确确定最优 k |
| 5 | 理论曲面与残差分析 | 测量值 vs 计算值散点图、残差直方图、理论曲面 |
| 6 | k 函数优化结果 | k 与 r、α 的关系，以及校正前后的角度对比 |
| 7 | θ 扰动优化结果 | 扰动项与 r、α 的关系，以及精度提升 |
| 8 | 校正函数优化对比 | 最终校正前后预测值与实测值的一致性 |

## 许可证

本项目采用 MIT 许可证。详见 LICENSE 文件。

## 联系方式

如有问题，请在 GitHub Issues 中提出。

（本README由DeepSeek、豆包辅助生成）



