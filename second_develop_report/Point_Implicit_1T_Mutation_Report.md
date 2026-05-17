# SU2 `1T + Mutation++` 点隐式化学源项报告

## 1. 报告目的

本文档面向当前仓库中的 `NEMO_EULER + Mutation++ + SINGLE_TEMPERATURE + POINT_IMPLICIT` 实现，系统说明以下内容：

- 单温度 `1T` Euler 框架下的控制方程；
- 与原始双温度 `2T` Euler 框架的根本差别；
- 有限体积半离散形式；
- 当前仓库中“流场显式、化学源项局部点隐式”的离散与线性化方式；
- 化学反应源项 Jacobian 的来源、分工与数学表达；
- 目前代码实现的适用范围、正确性结论与边界条件。

本文档的结论基于当前代码的静态审查。**本文不包含编译通过或算例收敛结果的运行验证**。

---

## 2. 适用范围

当前“真正的局部化学点隐式”只在如下条件下激活：

- 时间推进格式：`POINT_IMPLICIT`
- 主求解器：`NEMO_EULER`
- 热化学库：`MUTATIONPP`
- 温度模型：`SINGLE_TEMPERATURE`

对应代码：

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp:247`

即：

```cpp
return (config->GetKind_TimeIntScheme_Flow() == POINT_IMPLICIT) &&
       (config->GetKind_Solver() == MAIN_SOLVER::NEMO_EULER) &&
       (config->GetKind_FluidModel() == MUTATIONPP) &&
       (config->GetKind_NEMOTemperatureModel() == NEMO_TEMPERATURE_MODEL::SINGLE_TEMPERATURE);
```

这意味着：

- 本文讨论的是 **Euler 框架**；
- 本文讨论的是 **单温度 `1T`**；
- 本文讨论的是 **Mutation++ 路径**；
- 当前实现**不是** `NEMO_NS` 下的局部点隐式实现；
- 当前实现**不是** `2T` 下的局部点隐式实现。

---

## 3. `2T` 与 `1T` 的物理模型差别

### 3.1 `2T` 模型的本质

双温度 `2T` 非平衡模型的核心不是“有化学反应”，而是：

- 平动/转动模态使用一个温度 `T`；
- 振动/电子模态使用另一个温度 `Tve`；
- 振动/电子能量具有**独立演化方程**；
- 模态之间通过能量交换源项耦合。

因此，`2T` 的物理本质是：

1. 化学组成可能非平衡；
2. 热力学模态也可能非平衡；
3. `T` 与 `Tve` 不必相等。

### 3.2 `1T` 模型的本质

单温度 `1T` 模型并不是“没有振动能贡献”，而是：

- 不再给振动/电子能设置独立演化方程；
- 不再允许 `T` 与 `Tve` 独立偏离；
- 所有热模态统一由同一个温度 `T` 决定。

因此，`1T` 的本质是：

- `Tve = T`
- 振动能不再作为独立非平衡自由度演化
- 振动能、转动能、平动能仍可在热力学构成上区分，但它们都由同一个温度控制

当前代码中的直接体现：

- `1T` 时 `Tve` 被设置为 `T`
- `1T` 时第二套能量直接退化为第一套能量
- `1T` 时振动能交换源项为零

对应代码：

- `SU2_CFD/src/fluid/CMutationTCLib.cpp:156`
- `SU2_CFD/src/fluid/CMutationTCLib.cpp:165`
- `SU2_CFD/src/fluid/CMutationTCLib.cpp:232`

其中：

```cpp
if (single_temperature) {
  energies[1] = energies[0];
}
```

```cpp
SetTDStateRhosTTv(rhos, T, single_temperature ? T : val_T);
```

```cpp
if (single_temperature) {
  omega = 0.0;
  return omega;
}
```

---

## 4. Euler 框架下的控制方程

以下只讨论无粘 `Euler` 主体方程，不包含粘性通量。

### 4.1 `2T` 基础形式

`2T` 的守恒变量可写为：

\[
\mathbf{U}_{2T} =
\begin{bmatrix}
\rho_1,\dots,\rho_{N_s},
\rho u_1,\dots,\rho u_d,
\rho E,
\rho e_{ve}
\end{bmatrix}^{T}
\]

其中：

- `rho_s`：第 `s` 个组分的部分密度；
- `rho E`：总能量；
- `rho e_ve`：振动/电子能。

连续方程组可写为：

\[
\frac{\partial \rho_s}{\partial t} + \nabla \cdot (\rho_s \mathbf{u}) = \omega_s
\]

\[
\frac{\partial \rho \mathbf{u}}{\partial t} + \nabla \cdot (\rho \mathbf{u}\otimes\mathbf{u} + p\mathbf{I}) = 0
\]

\[
\frac{\partial \rho E}{\partial t} + \nabla \cdot \big[(\rho E + p)\mathbf{u}\big] = 0
\]

\[
\frac{\partial \rho e_{ve}}{\partial t} + \nabla \cdot (\rho e_{ve}\mathbf{u}) = \Omega_{ve}
\]

其中：

- `omega_s` 是化学反应引起的组分生成/消耗源项；
- `Omega_ve` 是平动/转动与振动/电子模态之间的能量交换源项。

### 4.2 当前 `1T` 模型的连续方程

在 `1T` 中，不再独立追踪 `rho e_ve` 的物理演化，因此物理上真正需要的是：

\[
\mathbf{U}_{1T} =
\begin{bmatrix}
\rho_1,\dots,\rho_{N_s},
\rho u_1,\dots,\rho u_d,
\rho E
\end{bmatrix}^{T}
\]

连续方程组为：

\[
\frac{\partial \rho_s}{\partial t} + \nabla \cdot (\rho_s \mathbf{u}) = \omega_s
\]

\[
\frac{\partial \rho \mathbf{u}}{\partial t} + \nabla \cdot (\rho \mathbf{u}\otimes\mathbf{u} + p\mathbf{I}) = 0
\]

\[
\frac{\partial \rho E}{\partial t} + \nabla \cdot \big[(\rho E + p)\mathbf{u}\big] = 0
\]

这里需要特别强调两点：

1. **单温度下仍然存在化学反应源项 `omega_s`**；
2. **单温度下总能量方程通常不再额外写一个独立“化学热源项”**，因为化学焓/内能变化已经通过组分变化和状态方程闭合进入总能量与温度关系中。

这也是当前代码中的实现方式：化学源项只进入组分方程行。

对应代码：

- `SU2_CFD/src/numerics/NEMO/NEMO_sources.cpp:97`

```cpp
for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++){
  residual[iSpecies] = ws[iSpecies] * Volume;
}
```

### 4.3 关于当前实现中的固定状态布局

当前 NEMO 代码框架为了复用原有结构，依然保留固定长度状态向量：

```cpp
U: [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
```

对应代码注释：

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp:101`

这意味着：

- 在实现层面，`1T` 并未把最后一个能量槽位从数据结构中彻底删除；
- 但在**物理闭合**上，`1T` 已经将 `Tve` 退化为 `T`，并关闭了独立的振动能交换源项；
- 因此当前 `1T` 的“物理有效模型”应理解为单温度热化学 Euler，而不是一个真正独立的双能量方程系统。

---

## 5. `1T` 热力学闭合

在 `1T` 中，总内能可形式化写为：

\[
e(T,\{Y_s\}) =
\sum_s Y_s
\Big(
e_{tr,s}(T) + e_{rot,s}(T) + e_{vib,s}(T) + e_{el,s}(T) + h^0_{f,s}
\Big)
\]

这里的关键不是是否把各模态概念上区分，而是：

\[
e_{tr,s},\ e_{rot,s},\ e_{vib,s},\ e_{el,s}
\quad \text{全部使用同一个温度 } T
\]

而在 `2T` 中，通常是：

\[
e_s = e_{tr,s}(T) + e_{rot,s}(T) + e_{vib,s}(T_{ve}) + e_{el,s}(T_{ve}) + h^0_{f,s}
\]

所以 `1T` 的本质不是“没有振动内能”，而是“振动内能不再独立非平衡演化”。

---

## 6. 有限体积半离散形式

设第 `i` 个控制体积体积为 `V_i`，则通用有限体积半离散形式可写为：

\[
V_i \frac{d \mathbf{U}_i}{dt}
+
\mathbf{R}_i^{conv}(\mathbf{U})
=
V_i \mathbf{S}_i(\mathbf{U}_i)
\]

对当前 `Euler + 1T` 问题，

- `R_i^{conv}` 为对流通量残差；
- `S_i` 为局部化学源项；
- 该源项只在组分方程上非零。

因此更紧凑地可写为：

\[
V_i \frac{d \mathbf{U}_i}{dt}
=
-\mathbf{R}_i(\mathbf{U})
+V_i \mathbf{S}_i(\mathbf{U}_i)
\]

其中 `R_i` 包含显式构造出的空间残差。

---

## 7. 当前点隐式方法的数值思想

### 7.1 基本思想

当前实现不是全局隐式，不组装整个流场的大稀疏 Jacobian，也不把对流通量 Jacobian 一起隐式化。

它做的是：

- 对流残差：显式处理；
- 化学源项：在每个网格点/单元处做局部线性化；
- 每个控制体独立求解一个小型稠密线性系统。

因此这是一种典型的：

**flow explicit + source local implicit**

即：

**流场显式、局部刚性源项隐式**

### 7.2 伪时间/时间推进形式

令 `Delta t_i` 表示当前单元使用的推进步长，则局部 Euler 线性化可写为：

\[
\frac{V_i}{\Delta t_i}
(\mathbf{U}_i^{n+1} - \mathbf{U}_i^{n})
=
-\mathbf{RHS}_i
+V_i \mathbf{S}_i(\mathbf{U}_i^{n+1})
\]

对源项在旧状态 `n` 处一阶线性化：

\[
\mathbf{S}_i(\mathbf{U}_i^{n+1})
\approx
\mathbf{S}_i(\mathbf{U}_i^{n})
+
\left(
\frac{\partial \mathbf{S}_i}{\partial \mathbf{U}}
\right)^n
\Delta \mathbf{U}_i
\]

其中：

\[
\Delta \mathbf{U}_i = \mathbf{U}_i^{n+1} - \mathbf{U}_i^{n}
\]

整理可得局部线性系统：

\[
\left(
\frac{V_i}{\Delta t_i}\mathbf{I}
-V_i \frac{\partial \mathbf{S}_i}{\partial \mathbf{U}}
\right)
\Delta \mathbf{U}_i
=
\mathbf{RHS}_i
\]

这正是当前代码的核心。

---

## 8. 当前代码中的局部线性系统

代码中局部线性系统直接写成：

\[
(\delta \mathbf{I} - \mathbf{J}_{chem}) \Delta \mathbf{U} = \mathbf{RHS}
\]

其中：

\[
\delta = \max(V_i/\Delta t_i,\ \varepsilon)
\]

对应代码：

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp:1008`

```cpp
/*--- Local dense point solves: (V/dt*I - J_chem) * dU = RHS. ---*/
```

以及：

```cpp
const su2double vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
const su2double delta = max(vol/dt, EPS);
```

```cpp
A[iVar][iVar] = delta;
...
A[iVar][jVar] -= Jchem[iVar][jVar];
```

因此当前局部点隐式矩阵就是：

\[
\mathbf{A}_i
=
\frac{V_i}{\Delta t_i}\mathbf{I}
-
\mathbf{J}_{chem,i}
\]

其中：

\[
\mathbf{J}_{chem,i}
=
\frac{\partial (V_i \mathbf{S}_i)}{\partial \mathbf{U}_i}
\]

### 8.1 残差右端

当前代码先使用显式残差和截断误差构造右端：

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp:981`

```cpp
LinSysRes[total_index] = -(LinSysRes[total_index] + local_Res_TruncError[iVar]);
```

因此右端本质上是：

\[
\mathbf{RHS}_i = -(\mathbf{R}_i + \mathbf{R}^{trunc}_i)
\]

### 8.2 求解与更新

局部系统用高斯消元求解：

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp:1045`

```cpp
Gauss_Elimination(A_rows, rhs, nVar);
```

得到 `Delta U` 后再乘欠松弛并更新：

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp:1095`

```cpp
nodes->AddSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint)*LinSysSol[iPoint*nVar+iVar]);
```

之后进行：

- 组分正性截断；
- 组分总和严格重标定。

对应代码：

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp:1106`

这一步是数值稳健性增强措施，不是控制方程本身的一部分。

---

## 9. 化学源项与 Jacobian 的结构

### 9.1 单元化学源项

当前 `1T` 化学源项向量可记为：

\[
\mathbf{S}(\mathbf{U}) =
\begin{bmatrix}
\omega_1 \\
\vdots \\
\omega_{N_s} \\
0 \\
\vdots \\
0
\end{bmatrix}
\]

其中：

- 前 `N_s` 行是组分生成/消耗率；
- 动量方程无化学源项；
- 总能量方程在当前 `1T` 保守总能量框架下不直接加入独立化学源项；
- 振动/电子能交换源项在 `1T` 下为零。

对应代码：

- `SU2_CFD/src/numerics/NEMO/NEMO_sources.cpp:97`
- `SU2_CFD/src/numerics/NEMO/NEMO_sources.cpp:135`

### 9.2 Jacobian 的总体链式法则

在 `1T` 下，化学反应速率依赖于：

- 组分密度 `rho_s`
- 单一温度 `T`

因此对守恒变量 `U` 的导数应写为：

\[
\frac{\partial \boldsymbol{\omega}}{\partial \mathbf{U}}
=
\frac{\partial \boldsymbol{\omega}}{\partial \boldsymbol{\rho}}\Big|_{T}
+
\frac{\partial \boldsymbol{\omega}}{\partial T}\Big|_{\rho}
\frac{\partial T}{\partial \mathbf{U}}
\]

这就是当前实现的核心分工：

- `Mutation++` 提供化学本征部分：
  - `omega`
  - `∂omega/∂rho|T`
  - `∂omega/∂T|rho`
- `SU2` 提供：
  - `dT/dU`

然后在 SU2 中完成最终链式组合。

---

## 10. `SU2` 侧的 `dT/dU`

当前单温度路径下，`dT/dU` 由 `CNEMOGas::ComputedTdU(...)` 给出：

- `SU2_CFD/src/fluid/CNEMOGas.cpp:253`

代码结构为：

- 对组分密度求导；
- 对动量分量求导；
- 对总能量求导；
- 对最后一个额外能量槽位在 `1T` 下给出零导数。

即：

\[
\frac{\partial T}{\partial (\rho E)} = \frac{1}{\rho C_{v,tr}}
\]

\[
\frac{\partial T}{\partial (\rho u_k)} = -\frac{u_k}{\rho C_{v,tr}}
\]

代码中还显式包含组分变化对温度的贡献：

\[
\frac{\partial T}{\partial \rho_s}
=
\frac{-e_{f,s} + \frac{1}{2}|\mathbf{u}|^2 + C_{v,tr,s}(T_{ref,s}-T)}
{\rho C_{v,tr}}
\]

从代码角度看，对应为：

```cpp
val_dTdU[iSpecies] = (...) / rhoCvtr;
val_dTdU[nSpecies+iDim] = -V[VEL_INDEX+iDim] / V[RHOCVTR_INDEX];
val_dTdU[nSpecies+nDim] = 1.0 / V[RHOCVTR_INDEX];
val_dTdU[nSpecies+nDim+1] = single_temperature ? 0.0 : (-1.0 / V[RHOCVTR_INDEX]);
```

这一步非常关键，因为：

- `Mutation++` 知道化学速率如何依赖 `T`；
- 但只有 `SU2` 知道自己的守恒变量 `U` 如何定义；
- 因此 `dT/dU` 必须在 `SU2` 侧完成。

---

## 11. `Mutation++` 提供的 Jacobian 组成

### 11.1 旧有能力：`jacobianRho`

`Mutation++` 原本已经能提供：

\[
\frac{\partial \boldsymbol{\omega}}{\partial \boldsymbol{\rho}}\Big|_T
\]

对应接口：

- `mix->jacobianRho(...)`

这条路径在非 `1T` 真点隐式分支或其他隐式路径中仍然存在。

### 11.2 新增能力：`netProductionRatesJacobian1T`

当前为了 `1T` 真点隐式路径，新增了接口：

- `subprojects/Mutationpp/src/kinetics/Kinetics.h`
- `subprojects/Mutationpp/src/kinetics/Kinetics.cpp:370`

接口功能是一次性返回：

- `wdot`
- `∂omega/∂rho|T`
- `∂omega/∂T|rho`

即：

```cpp
void netProductionRatesJacobian1T(double* p_wdot, double* p_jac, double* p_dwdt);
```

---

## 12. `Mutation++` 中 `∂omega/∂rho|T` 的来源

在 `netProductionRatesJacobian1T(...)` 中，先构造：

- 正反应速率系数；
- 正反应/逆反应反应进度；
- 第三体修正；
- 物种浓度。

然后调用：

```cpp
m_jacobian.computeJacobian(mp_ropf, mp_ropb, conc.data(), p_jac);
```

对应代码：

- `subprojects/Mutationpp/src/kinetics/Kinetics.cpp:412`

这一步给出的就是：

\[
\frac{\partial \boldsymbol{\omega}}{\partial \boldsymbol{\rho}}\Big|_T
\]

它是 `Mutation++` 内部原生 Jacobian 管理器计算出来的，不是 SU2 自己再推导的。

---

## 13. `Mutation++` 中 `∂omega/∂T|rho` 的来源

### 13.1 基本公式

对第 `r` 个反应，记净反应进度为：

\[
\mathcal{R}_r = \mathcal{R}_{f,r} - \mathcal{R}_{b,r}
\]

则对温度求导：

\[
\frac{\partial \mathcal{R}_r}{\partial T}
=
\frac{\partial \mathcal{R}_{f,r}}{\partial T}
-
\frac{\partial \mathcal{R}_{b,r}}{\partial T}
\]

若把浓度固定，则有：

\[
\frac{\partial \mathcal{R}_{f,r}}{\partial T}
=
\frac{\partial \ln k_{f,r}}{\partial T}\mathcal{R}_{f,r}
\]

\[
\frac{\partial \mathcal{R}_{b,r}}{\partial T}
=
\frac{\partial \ln k_{b,r}}{\partial T}\mathcal{R}_{b,r}
\]

因此：

\[
\frac{\partial \mathcal{R}_r}{\partial T}
=
\frac{\partial \ln k_{f,r}}{\partial T}\mathcal{R}_{f,r}
-
\frac{\partial \ln k_{b,r}}{\partial T}\mathcal{R}_{b,r}
\]

### 13.2 Arrhenius 正反应项

当前实现假定反应率型为 Arrhenius。代码中得到：

\[
\frac{\partial \ln k_f}{\partial T}
=
\frac{n + T_c/T}{T}
\]

对应代码：

- `subprojects/Mutationpp/src/kinetics/Kinetics.cpp:444`

```cpp
const double dlnkf_dT = (p_rate->n() + p_rate->T()/T) / T;
```

### 13.3 逆反应项

对于可逆反应：

\[
\ln k_b = \ln k_f - \ln K_{eq}
\]

因此：

\[
\frac{\partial \ln k_b}{\partial T}
=
\frac{\partial \ln k_f}{\partial T}
-
\frac{\partial \ln K_{eq}}{\partial T}
\]

当前实现使用：

\[
\frac{\partial \ln K_{eq}}{\partial T}
=
\frac{\Delta(h^0/RT)-\Delta \nu}{T}
\]

其中：

- `Delta(h^0/RT)` 来自产物与反应物的 `speciesHOverRT` 差；
- `Delta nu` 为化学计量数总和差。

对应代码：

- `subprojects/Mutationpp/src/kinetics/Kinetics.cpp:449`

```cpp
const double dlnkeq_dT = (delta_h_over_rt - delta_nu) / T;
dlnkb_dT = dlnkf_dT - dlnkeq_dT;
```

### 13.4 反应进度对温度导数

因此代码中构造：

\[
\frac{\partial \mathcal{R}_r}{\partial T}
=
\frac{\partial \ln k_{f,r}}{\partial T}\mathcal{R}_{f,r}
-
\frac{\partial \ln k_{b,r}}{\partial T}\mathcal{R}_{b,r}
\]

对应代码：

- `subprojects/Mutationpp/src/kinetics/Kinetics.cpp:466`

```cpp
const double drop_dT =
    dlnkf_dT * mp_ropf[i_rxn] - dlnkb_dT * mp_ropb[i_rxn];
```

### 13.5 从反应进度导数累加到物种生成率导数

若 `nu_sr` 是第 `r` 个反应对第 `s` 个组分的化学计量贡献，则：

\[
\frac{\partial \omega_s}{\partial T}
=
M_s \sum_r \nu_{sr}\frac{\partial \mathcal{R}_r}{\partial T}
\]

当前实现正是通过对反应物减、对产物加的方式逐反应累加：

- `subprojects/Mutationpp/src/kinetics/Kinetics.cpp:472`

这是解析导数，不再是旧版本那种有限差分温度扰动近似。

### 13.6 第三体项的温度导数

当前实现中第三体因子先乘到正反应/逆反应进度上：

- `subprojects/Mutationpp/src/kinetics/Kinetics.cpp:409`

在当前 `1T`、固定组分密度条件下，第三体因子的直接温度导数不单独显式加入，这是合理的，因为第三体项本身依赖的是浓度组合，而不是显式的温度函数。

---

## 14. SU2 中最终 Jacobian 的装配

当前 `1T` 单温度点隐式路径在 `CMutationTCLib::ComputeNetProductionRates(...)` 中执行最终装配：

- `SU2_CFD/src/fluid/CMutationTCLib.cpp:178`

装配过程分两步：

### 14.1 先装配物种密度块

```cpp
val_jacobian[iSpecies][jSpecies] += jac_species[iSpecies*nSpecies + jSpecies];
```

这对应：

\[
\frac{\partial \omega_i}{\partial \rho_j}\Big|_T
\]

### 14.2 再加链式温度项

```cpp
val_jacobian[iSpecies][jVar] += dwdT[iSpecies] * dTdU[jVar];
```

这对应：

\[
\frac{\partial \omega_i}{\partial U_j}
\leftarrow
\frac{\partial \omega_i}{\partial U_j}
+
\frac{\partial \omega_i}{\partial T}
\frac{\partial T}{\partial U_j}
\]

因此最终 `SU2` 使用的 Jacobian 是：

\[
\boxed{
\frac{\partial \boldsymbol{\omega}}{\partial \mathbf{U}}
=
\frac{\partial \boldsymbol{\omega}}{\partial \boldsymbol{\rho}}\Big|_{T}
+
\frac{\partial \boldsymbol{\omega}}{\partial T}\Big|_{\rho}
\frac{\partial T}{\partial \mathbf{U}}
}
\]

这是当前实现中最核心的数学关系。

---

## 15. 单元体积缩放

`Mutation++` 返回的是单位体积意义下的源项导数，而单元离散需要的是：

\[
\frac{\partial (V_i \mathbf{S}_i)}{\partial \mathbf{U}_i}
\]

因此在 `AssembleLocalChemJacobian(...)` 中，SU2 会再乘单元体积：

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp:279`

```cpp
local_jac[iVar][jVar] *= volume;
```

于是局部系统中的 `Jchem` 实际上是：

\[
\mathbf{J}_{chem,i}
=
V_i
\frac{\partial \mathbf{S}_i}{\partial \mathbf{U}_i}
\]

---

## 16. 当前点隐式的作用

当前局部点隐式的作用不是把整个流场都做成全隐式，而是专门抑制**化学刚性**。

其主要收益是：

- 改善刚性化学源项导致的数值稳定性问题；
- 减小显式化学推进对步长的苛刻限制；
- 在定常伪时间推进中允许更积极的 CFL 或伪时间步；
- 在非定常内迭代中改善每个物理时间步内部的收敛性。

但它**不是**：

- 对物理时间步无限制放大的许可；
- 对流通量全隐式；
- 全耦合 Newton-Krylov 全局非线性隐式。

所以应当把它理解为：

**局部源项稳定化 / 局部刚性预处理**

而不是完整全系统隐式化。

---

## 17. 对定常与非定常的适用性

从代码入口看，`POINT_IMPLICIT` 是统一调用的，不限定只能用于定常：

- `SU2_CFD/src/integration/CIntegration.cpp`

因此从实现上讲：

- 可以用于定常伪时间推进；
- 也可以用于非定常推进中的内层迭代。

但从数值意义上讲应区分：

### 17.1 定常问题

对定常问题，局部点隐式通常可以：

- 提高稳定性；
- 加快收敛；
- 允许更大的伪时间步。

### 17.2 非定常问题

对非定常问题，局部点隐式可以：

- 改善内层求解稳定性；
- 改善刚性源项处理；

但不能改变这样一个事实：

- **物理时间步长仍然受时间精度要求约束。**

因此：

- 它可以提高稳定性；
- 但不能把真实非定常时间积分的精度约束“消掉”。

---

## 18. 当前实现的正确性结论

在当前限定范围内，即：

- `Mutation++`
- `1T`
- 非电离
- Arrhenius 反应率
- `NEMO_EULER`
- `POINT_IMPLICIT`

当前实现从静态审查角度看是**物理自洽、数值结构正确、代码路径闭合**的。

主要依据如下：

1. `1T` 下化学源项只作用于组分方程，这与保守总能量框架一致；
2. `1T` 下不再独立求解振动能交换源项，符合 `Tve = T` 的模型退化；
3. 局部点隐式矩阵形式为标准的
   \[
   (V/\Delta t I - J_{chem})\Delta U = RHS
   \]
4. `∂omega/∂U` 的链式分解正确；
5. `∂omega/∂rho|T` 由 `Mutation++` 内部 Jacobian 管理器给出；
6. `∂omega/∂T|rho` 已由解析公式替代旧有限差分近似；
7. 最终 Jacobian 装配与体积缩放关系正确。

---

## 19. 当前实现的边界与限制

### 19.1 仅覆盖 `1T`

新解析 `dwdT` 路径只针对：

- `nEnergyEqns() == 1`

对应代码：

- `subprojects/Mutationpp/src/kinetics/Kinetics.cpp:376`

### 19.2 仅覆盖 Arrhenius 速率

当前 `Mutation++` 的解析 `∂omega/∂T` 实现要求：

- 反应率律可转型为 `Arrhenius`

否则直接抛出 `NotImplementedError`。

### 19.3 当前真点隐式只对 `NEMO_EULER` 生效

`NEMO_NS` 目前并未切入这条“真正的局部化学点隐式”路径。

### 19.4 未包含本地编译/算例验证

本文档结论来自代码静态审查，不等价于：

- 已完成编译通过验证；
- 已完成单元测试；
- 已完成基准算例对照。

---

## 20. 建议的后续验证工作

为把静态正确性进一步提升到工程可信度，建议后续做以下验证：

1. 编译级验证  
   确认 `Mutation++` 新接口、SU2 调用链、链接过程全部通过。

2. Jacobian 一致性验证  
   在单点状态上比较：
   - 解析 `∂omega/∂T`
   - 有限差分 `∂omega/∂T`
   的误差。

3. 局部线性化验证  
   对比：
   - `J * deltaU`
   - `omega(U+deltaU)-omega(U)`
   的一致性。

4. 算例验证  
   对 `1T` 非电离反应流算例比较：
   - 显式化学
   - 点隐式化学
   的稳定性、收敛性和结果一致性。

5. 边界状态验证  
   检查高温、低密度、稀薄激波后状态下 Jacobian 与更新的稳健性。

---

## 21. 总结

当前仓库中的 `1T + Mutation++ + POINT_IMPLICIT` 路径，本质上是：

- 在单温度热化学 Euler 模型下；
- 用显式方式处理流场对流残差；
- 用局部点隐式方式处理刚性化学源项；
- 其核心线性系统为
  \[
  \left(\frac{V_i}{\Delta t_i}I - J_{chem,i}\right)\Delta U_i = RHS_i
  \]
- 其中
  \[
  J_{chem}
  =
  V_i
  \left[
  \frac{\partial \omega}{\partial \rho}\Big|_T
  +
  \frac{\partial \omega}{\partial T}\Big|_\rho
  \frac{\partial T}{\partial U}
  \right]
  \]

当前这套实现的重要进展在于：

- `∂omega/∂T` 已不再依赖有限差分，而是切换为 `Mutation++` 内部解析构造；
- `SU2` 与 `Mutation++` 的职责边界清晰：
  - `Mutation++` 负责化学本征导数；
  - `SU2` 负责守恒变量到温度的链式映射；
- 在当前限定范围内，这是一套物理和数值上都合理的 `1T` 局部点隐式实现。

