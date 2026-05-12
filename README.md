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

# 算法步骤
**1. 常数 k 拟合**
基于理论解析式，全局搜索最优常数 $k$：
$$
\theta_{\rm theory} = \mathrm{atan2}\left( 1 \pm \sqrt{1 - k r \left(2\sin\alpha + k r \cos^2\alpha\right)},\; k r \cos\alpha \right) \times \frac{180}{\pi}
$$

**2. k 函数优化（步骤6）**
将 $k$ 扩展为关于 $r,\alpha$ 的二次耦合函数，拟合全局参数：
$$
k = k_{\rm base} \cdot \left(1 + a r + b \alpha + c r\alpha + d r^2 + e \alpha^2\right)
$$
- 输出参数：$k_{\rm base}、a、b、c、d、e$
- 保存文件：`step6_results.mat`

**3. θ 扰动优化（步骤7）**
引入多物理项扰动，修正理论角度固有偏差：
$$
\theta_{\rm pert} = \theta_{\rm theory} + A_1 + A_2 r + A_3 \sin\alpha + A_4 \sin2\alpha + A_5 r\sin\alpha
$$
- 输出参数：A1~A5
- 保存文件：`step7_results.mat`

**4. 分段校正（步骤8）**
以 $\boldsymbol{45^\circ}$ 为阈值，采用两组线性模型分段校正角度误差 $\Delta\theta$
- 输出两套校正系数
- 保存文件：`step8_results.mat`

---

# 最终预测公式
将上述步骤优化得到的参数代入，即可计算任意 $(r,\alpha)$ 对应的预测角度 $\theta_{\rm pred}$。

## 1. 原始理论角度解析式
$$
\theta_{\rm theory} = \mathrm{atan2}\left( 1 \pm \sqrt{1 - k r \left(2\sin\alpha + k r \cos^2\alpha\right)},\; k r \cos\alpha \right) \times \frac{180}{\pi}
$$

## 2. k 参数化表达式（取自 `step6_results.mat`）
$$
k = k_{\rm base} \cdot \left(1 + a r + b \alpha + c r\alpha + d r^2 + e \alpha^2\right)
$$

## 3. 扰动后角度表达式（取自 `step7_results.mat`）
$$
\theta_{\rm pert} = \theta_{\rm theory} + A_1 + A_2 r + A_3 \sin\alpha + A_4 \sin2\alpha + A_5 r\sin\alpha
$$

## 4. 分段校正项表达式（取自 `step8_results.mat`）
- 情况1：$\boldsymbol{\theta_{\rm pert} \le 45^\circ}$
$$
\begin{aligned}
\Delta\theta &= c_0 + c_1\theta_{\rm pert} + c_2\theta_{\rm pert}^2 + c_3\left(\frac{r}{300}\right) + c_4\left(\frac{r}{300}\right)^2 \\
&\quad + c_5\sin\alpha + c_6\cos\alpha + c_7\sin2\alpha + c_8\left(\frac{r}{300}\right)\sin\alpha + c_9\left(\frac{\alpha}{90}\right)
\end{aligned}
$$

- 情况2：$\boldsymbol{\theta_{\rm pert} > 45^\circ}$
$$
\Delta\theta = d_0 + d_1\theta_{\rm pert} + d_2\theta_{\rm pert}^2 + d_3\left(\frac{r}{300}\right) + d_4\left(\frac{r}{300}\right)^2 + d_5\sin\alpha + d_6\sin2\alpha
$$

## 最终预测值
$$
\boldsymbol{\theta_{\rm pred} = \theta_{\rm pert} + \Delta\theta}

## 输出图形说明

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



