# SU2 NEMO 二次开发报告

## 1. 报告目的

本文档面向当前仓库后续的 NEMO 框架二次开发工作，目标不是做一般性的 SU2 介绍，而是建立一份可直接指导后续开发的技术地图，重点说明：

- 当前分支的整体代码框架与运行链路。
- NEMO 模块的控制方程、物理闭合、数值方法与代码落点。
- 当前分支中已确认存在的二次开发内容及其真实实现范围。
- 后续继续扩展物理模型、化学模型、数值方法和网格自适应时，应该从哪里切入，以及哪些位置最脆弱。

本文档基于当前仓库代码直接阅读得到，不依赖 git 历史。

---

## 2. 仓库整体框架

### 2.1 宏观执行链路

当前分支仍然保持 SU2 的标准主链路：

`config -> driver -> iteration -> solver -> numerics/fluid/variables -> output`

关键入口如下：

- 主程序入口：[SU2_CFD/src/SU2_CFD.cpp](SU2_CFD/src/SU2_CFD.cpp)
- 求解器工厂：[SU2_CFD/src/solvers/CSolverFactory.cpp](SU2_CFD/src/solvers/CSolverFactory.cpp)
- 主迭代控制：[SU2_CFD/src/iteration/CFluidIteration.cpp](SU2_CFD/src/iteration/CFluidIteration.cpp)

其逻辑为：

1. 读取配置文件并建立 `CConfig`。
2. 根据求解器类型选择 driver。
3. driver 创建几何、变量、求解器、离散算子和输出对象。
4. 主迭代层调用流场求解器。
5. 流场求解器进一步组织：
   - 原始变量恢复
   - 对流残差
   - 黏性残差
   - 源项残差
   - 时间推进
   - 输出与收敛监控

### 2.2 NEMO 在整体框架中的位置

NEMO 不是通过“可压缩流 + 若干附加标量”的方式临时拼接，而是直接作为主求解器接入框架。

在 [SU2_CFD/src/solvers/CSolverFactory.cpp](SU2_CFD/src/solvers/CSolverFactory.cpp) 中已确认：

- `MAIN_SOLVER::NEMO_EULER -> CNEMOEulerSolver`
- `MAIN_SOLVER::NEMO_NAVIER_STOKES -> CNEMONSSolver`
- `MAIN_SOLVER::NEMO_RANS -> CNEMONSSolver` 并叠加湍流/转捩模块

因此，NEMO 的方程组、状态恢复、源项、边界条件、输出和网格自适应接口都已经进入 SU2 的主干求解链条。

---

## 3. NEMO 的物理模型与状态设计

### 3.1 NEMO 的状态变量

在 [SU2_CFD/src/variables/CNEMOEulerVariable.cpp](SU2_CFD/src/variables/CNEMOEulerVariable.cpp) 和 [SU2_CFD/src/solvers/CNEMOEulerSolver.cpp](SU2_CFD/src/solvers/CNEMOEulerSolver.cpp) 中可确认，NEMO 的守恒量布局为：

`U = [rho_1, ..., rho_Ns, rho*u, rho*v, rho*w, rho*E, rho*Eve]^T`

其中：

- `rho_s`：第 `s` 个组分的偏密度
- `rho*u_i`：动量
- `rho*E`：总能
- `rho*Eve`：振动-电子能

原始变量布局为：

`V = [rho_1, ..., rho_Ns, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T`

这说明当前 NEMO 框架描述的是：

- 多组分气体
- 热化学非平衡流
- 平动/转动温度 `T`
- 振动/电子温度 `Tve`

如果是单温模型，则 `Tve = T`，但代码结构仍保留 `rhoEve` 方程槽位和相关输出接口。

### 3.2 NEMO 的热力学闭合接口

NEMO 的热力学和化学闭合由抽象基类 [SU2_CFD/include/fluid/CNEMOGas.hpp](SU2_CFD/include/fluid/CNEMOGas.hpp) 统一定义，主要负责：

- 根据 `rho_s, T, Tve` 建立热力学状态
- 计算压强、声速、比热比
- 计算组分焓、混合物能量
- 计算化学生成率
- 计算振动弛豫源项
- 计算输运系数
- 提供 `dPdU`、`dTdU`、`dTvedU`

具体后端有两套：

- Mutation++ 接口：[SU2_CFD/src/fluid/CMutationTCLib.cpp](SU2_CFD/src/fluid/CMutationTCLib.cpp)
- SU2 自带 nonequilibrium 化学库：[SU2_CFD/src/fluid/CSU2TCLib.cpp](SU2_CFD/src/fluid/CSU2TCLib.cpp)

这意味着 NEMO 求解器本身并不直接固化具体气体模型，而是通过统一接口调用不同热化学后端。

---

## 4. NEMO 的数值主链路

### 4.1 守恒量到原始变量的恢复

`Cons2PrimVar` 是 NEMO 中最关键的数值节点，位置在：

- [SU2_CFD/src/variables/CNEMOEulerVariable.cpp](SU2_CFD/src/variables/CNEMOEulerVariable.cpp)

其核心工作包括：

1. 对各组分偏密度进行正性保护。
2. 累加得到总密度。
3. 由守恒量恢复速度。
4. 根据 `rhoE`、`rhoEve`、动能和组分密度，通过热力学库反解温度。
5. 更新：
   - `P`
   - `a`
   - `h`
   - `rhoCvtr`
   - `rhoCvve`
   - `eves`
   - `Cvves`
   - `dPdU`
   - `dTdU`
   - `dTvedU`

这一步的物理意义是：离散系统实际推进的是守恒量，但通量、源项、边界条件和大多数物性闭合都需要原始变量和热力学导数，因此每一步都要通过状态反演把数值状态重新映射回物理状态。

### 4.2 对流项离散

NEMO 的对流离散位于：

- [SU2_CFD/src/numerics/NEMO/CNEMONumerics.cpp](SU2_CFD/src/numerics/NEMO/CNEMONumerics.cpp)
- [SU2_CFD/src/numerics/NEMO/convection/roe.cpp](SU2_CFD/src/numerics/NEMO/convection/roe.cpp)
- [SU2_CFD/src/numerics/NEMO/convection/ausm_slau.cpp](SU2_CFD/src/numerics/NEMO/convection/ausm_slau.cpp)
- [SU2_CFD/src/numerics/NEMO/convection/msw.cpp](SU2_CFD/src/numerics/NEMO/convection/msw.cpp)

实现特征如下：

- Roe 方法使用完整特征分解。
- AUSM/SLAU 类型方法使用质量流量分裂与压强分裂。
- MSW 使用加权星状态来增强间断附近鲁棒性。

重要的是，NEMO 的这些格式并不是“主流体 + 被动组分输运”的简单叠加，而是把以下信息真正纳入了 Jacobian 和特征结构：

- 组分偏密度变化
- 压强对守恒量的导数
- 声速变化
- 振动能耦合
- 单温与双温模型差异

这是必要的，因为高温化学非平衡流中，激波后状态、组分变化和能量模态之间是强耦合的，若只把组分当被动标量，通量线性化会明显失真。

### 4.3 源项离散

NEMO 源项位于：

- [SU2_CFD/src/numerics/NEMO/NEMO_sources.cpp](SU2_CFD/src/numerics/NEMO/NEMO_sources.cpp)

分三类：

- 化学反应源项 `ComputeChemistry`
- 振动弛豫源项 `ComputeVibRelaxation`
- 轴对称源项 `ComputeAxisymmetric`

在 [SU2_CFD/src/solvers/CNEMOEulerSolver.cpp](SU2_CFD/src/solvers/CNEMOEulerSolver.cpp) 的 `Source_Residual` 中，源项被逐点加入线性系统，并在隐式模式下同步写入对角 Jacobian 块。

---

## 5. NEMO 的物理规律、数学描述与代码实现

### 5.1 有限速率化学

SU2 自带 nonequilibrium 化学实现位于：

- [SU2_CFD/src/fluid/CSU2TCLib.cpp](SU2_CFD/src/fluid/CSU2TCLib.cpp)

核心函数：

- `ComputeNetProductionRates`
- `ChemistryJacobian`
- `ComputeKeqConstants`

数学思想如下：

1. 对每条反应定义控制反应速率的温度。
2. 正反应率按 Arrhenius 形式计算。
3. 逆反应率通过平衡常数 `Keq` 得到。
4. 根据质量作用定律计算正向和逆向反应进度。
5. 按化学计量关系汇总到各组分生成率 `ws`。

当前实现中还有两个重要数值处理：

- 反应温度采用平滑修正温度 `Thf / Thb`
- 平衡常数通过密度相关表插值获得系数

平滑温度的目的不是改变物理模型，而是缓和源项刚性，避免在低温或过渡区出现过强的局部刚性，提升收敛性。

### 5.2 化学 Jacobian

`ChemistryJacobian` 在 [SU2_CFD/src/fluid/CSU2TCLib.cpp](SU2_CFD/src/fluid/CSU2TCLib.cpp) 中将化学源项对状态的导数显式展开，主要包括：

- `dkf/dU`
- `dkb/dU`
- `dRf/dU`
- `dRb/dU`
- 振动能源项中的组分能贡献导数

这一步的意义在于：

- 若采用隐式推进，化学刚性必须通过 Jacobian 才能稳定地进入离散系统。
- 若不考虑这些导数，强放热、强解离或强复合反应会导致显式或弱隐式方法稳定性急剧恶化。

### 5.3 振动弛豫

振动-平动能量交换在 `ComputeEveSourceTerm` 中实现，位于：

- [SU2_CFD/src/fluid/CSU2TCLib.cpp](SU2_CFD/src/fluid/CSU2TCLib.cpp)

所用模型可概括为：

- Landau-Teller 形式
- Millikan-White 弛豫时间
- Park 限幅碰撞截面

并考虑了：

- 振动温度向平动温度弛豫
- 化学反应引起的振动能变化

这在高温空气动力学中是合理且常见的工程模型，因为：

- 平动/转动模态响应快
- 振动模态响应相对慢
- 激波层和高焓流中，振动非平衡对温度场、声速和反应速率都会产生重要反馈

### 5.4 边界条件

NEMO 边界条件主实现位于：

- [SU2_CFD/src/solvers/CNEMOEulerSolver.cpp](SU2_CFD/src/solvers/CNEMOEulerSolver.cpp)
- [SU2_CFD/src/solvers/CNEMONSSolver.cpp](SU2_CFD/src/solvers/CNEMONSSolver.cpp)

已确认情况：

- `BC_Far_Field` 可用
- `BC_Supersonic_Inlet` 可用
- `BC_Supersonic_Outlet` 可用
- `BC_Outlet` 含特征外推逻辑
- `BC_Inlet` 明确标记为 `Not operational in NEMO`

壁面条件方面：

- 非催化热流壁可用
- 非催化等温壁可用
- 催化等温壁至少部分实现
- 催化热流壁和 wall function 在 NEMO 中并不成熟

结论是：当前 NEMO 的主体求解链路成熟，但壁面催化和部分入口条件仍然是脆弱区域，后续开发时不能默认其完备性。

---

## 6. 已确认存在的二次开发内容

本节只陈述已经在代码中确认存在的事实。

### 6.1 网格自适应分支已并入当前分支

已确认存在 Python 层的 MMG 外循环自适应流程：

- [SU2_PY/mesh_adaptation.py](SU2_PY/mesh_adaptation.py)
- [SU2_PY/SU2/adap/mmg.py](SU2_PY/SU2/adap/mmg.py)
- [SU2_PY/SU2/adap/tools.py](SU2_PY/SU2/adap/tools.py)

其工作模式不是求解器内部动态重构网格，而是：

1. 运行 CFD。
2. 输出场变量与 metric。
3. 调用 MMG 进行外部重构网格。
4. 插值或重启后继续求解。

NEMO 侧的适配接口在：

- [SU2_CFD/include/solvers/CNEMOEulerSolver.hpp](SU2_CFD/include/solvers/CNEMOEulerSolver.hpp)

已实现的 NEMO 自适应传感量包括：

- `MACH`
- `PRESSURE`
- `TEMPERATURE`
- `TEMPERATURE_VE`
- `ENERGY`
- `ENERGY_VE`

说明你们嵌入的自适应功能已经接入 NEMO 变量输出链路，而不是仅存在于普通可压缩流模块中。

### 6.2 NEMO 单温模型已系统实现

已确认配置入口：

- [Common/include/option_structure.hpp](Common/include/option_structure.hpp)
- [Common/src/CConfig.cpp](Common/src/CConfig.cpp)

已确认实现贯穿以下层次：

- 配置层：`NEMO_TEMPERATURE_MODEL = SINGLE_TEMPERATURE`
- 物性层：Mutation++ 选择 `ChemNonEq1T`
- 状态恢复层：`Tve = T`，并重构一致的 `rhoEve`
- 数值层：单温模式下屏蔽或合并相关振动项
- 源项层：单温模式下振动弛豫项为零
- 输出层：`TEMPERATURE_VE`、`ENERGY_VE` 对应单温逻辑

这说明单温模型不是表面开关，而是完整贯穿：

`config -> fluid -> variables -> numerics -> sources -> output`

### 6.3 Mutation++ 化学源项局部点隐式方法已实现

核心位置：

- [SU2_CFD/src/solvers/CNEMOEulerSolver.cpp](SU2_CFD/src/solvers/CNEMOEulerSolver.cpp)
- [SU2_CFD/src/fluid/CMutationTCLib.cpp](SU2_CFD/src/fluid/CMutationTCLib.cpp)

已确认的实现逻辑如下：

1. 只有在以下条件同时满足时，进入“真点隐式模式”：
   - `POINT_IMPLICIT`
   - `NEMO_EULER`
   - `MUTATIONPP`
   - `SINGLE_TEMPERATURE`

2. 此模式下不使用普通的全局隐式装配，而是逐点建立局部稠密块：

   `(V/dt * I - Jchem) * dU = RHS`

3. `Jchem` 来自局部化学 Jacobian 装配：
   - Mutation++ 的 `jacobianRho`
   - 外加有限差分获得的 `domega/dT * dT/dU` 链式项

4. 更新完成后执行：
   - species 正性截断
   - 组分密度严格重标定
   - 欠松弛因子限制

这说明你们实现的不是普通意义上的“点迭代器”，而是针对刚性化学源项的局部隐式稳定化机制。

其数值动机很明确：

- 高温有限速率化学往往局部极刚。
- 若上全局强耦合隐式，代价过高。
- 采用局部块隐式可以以较低成本显著增强稳定性。

### 6.4 SU2_NONEQ 自定义化学机理已扩展

配置入口已确认：

- `SU2_NONEQ_CUSTOM_CHEMISTRY`
- `SU2_NONEQ_SPECIES_TABLE`
- `SU2_NONEQ_ELECTRONIC_STATES_TABLE`
- `SU2_NONEQ_REACTIONS_TABLE`
- `SU2_NONEQ_TRANSPORT_TABLE`
- `SU2_NONEQ_KEQ_TABLE`

对应代码位于：

- [Common/src/CConfig.cpp](Common/src/CConfig.cpp)
- [SU2_CFD/src/fluid/CSU2TCLib.cpp](SU2_CFD/src/fluid/CSU2TCLib.cpp)

这表明你们并非只修改了少量反应参数，而是把 SU2 自带 nonequilibrium 化学库扩展成了可由外部表驱动的框架，能够加载：

- 组分数据
- 电子能级
- 反应机理
- 平衡常数系数
- 输运性质

这是一项结构性改动，后续任何新的反应机理扩展都可以继续沿这条路径推进。

---

## 7. NEMO 代码的后续开发入口图

### 7.1 若要修改控制方程或状态变量定义

首要入口：

- [SU2_CFD/src/variables/CNEMOEulerVariable.cpp](SU2_CFD/src/variables/CNEMOEulerVariable.cpp)
- [SU2_CFD/include/fluid/CNEMOGas.hpp](SU2_CFD/include/fluid/CNEMOGas.hpp)
- [SU2_CFD/src/solvers/CNEMOEulerSolver.cpp](SU2_CFD/src/solvers/CNEMOEulerSolver.cpp)

需要同步检查：

- 原始变量索引
- `Cons2PrimVar`
- `RecomputeConservativeVector`
- `CheckNonPhys`
- 输出字段映射

原因是 NEMO 的很多模块都假定当前状态布局固定，一旦状态变量改变，最先出问题的通常不是主残差，而是：

- 变量恢复
- Jacobian 导数
- 输出/重启文件兼容性

### 7.2 若要修改化学机理

Mutation++ 路径：

- [SU2_CFD/src/fluid/CMutationTCLib.cpp](SU2_CFD/src/fluid/CMutationTCLib.cpp)

SU2_NONEQ 路径：

- [SU2_CFD/src/fluid/CSU2TCLib.cpp](SU2_CFD/src/fluid/CSU2TCLib.cpp)

若仅新增机理数据，优先走表驱动路径。

若要修改公式级实现，应重点关注：

- `ComputeNetProductionRates`
- `ChemistryJacobian`
- `ComputeKeqConstants`
- `ComputeEveSourceTerm`
- `GetEveSourceTermJacobian`

### 7.3 若要修改刚性源项离散或时间推进

重点位置：

- [SU2_CFD/src/solvers/CNEMOEulerSolver.cpp](SU2_CFD/src/solvers/CNEMOEulerSolver.cpp)
- [SU2_CFD/src/numerics/NEMO/NEMO_sources.cpp](SU2_CFD/src/numerics/NEMO/NEMO_sources.cpp)

特别需要关注：

- `Source_Residual`
- `PointImplicit_Iteration`
- `AssembleLocalChemJacobian`
- `CompleteImplicitIteration`
- 欠松弛与正性保护

这部分直接决定刚性问题是否稳定，是后续数值开发最敏感的区域之一。

### 7.4 若要修改边界条件

重点位置：

- [SU2_CFD/src/solvers/CNEMOEulerSolver.cpp](SU2_CFD/src/solvers/CNEMOEulerSolver.cpp)
- [SU2_CFD/src/solvers/CNEMONSSolver.cpp](SU2_CFD/src/solvers/CNEMONSSolver.cpp)

注意事项：

- 入口边界在 NEMO 中并不完全统一成熟。
- 催化壁相关功能成熟度不均。
- 修改边界条件时必须同步检查热力学状态重构是否一致。

### 7.5 若要修改网格自适应指标

重点位置：

- [SU2_CFD/include/solvers/CNEMOEulerSolver.hpp](SU2_CFD/include/solvers/CNEMOEulerSolver.hpp)
- [SU2_CFD/src/drivers/CSinglezoneDriver.cpp](SU2_CFD/src/drivers/CSinglezoneDriver.cpp)
- [SU2_PY/SU2/adap/mmg.py](SU2_PY/SU2/adap/mmg.py)
- [SU2_PY/SU2/adap/tools.py](SU2_PY/SU2/adap/tools.py)

流程是：

1. C++ 侧定义 NEMO 要输出的 adaptation auxiliary variables。
2. C++ 侧生成梯度和 metric。
3. Python 侧读取输出文件。
4. MMG 使用 metric 进行网格调整。

如果后续要引入更适合高温流的传感量，比如反应进度变量、组分梯度范数、弛豫非平衡度等，推荐从这里切入。

---

## 8. 当前分支的成熟区域与高风险区域

### 8.1 相对成熟的区域

- NEMO 主求解链路
- 单温模型的全链路贯通
- Mutation++ 单温点隐式化学稳定化
- SU2_NONEQ 自定义化学表驱动扩展
- MMG 外循环自适应与 NEMO 变量输出耦合

### 8.2 高风险区域

- NEMO 壁面催化边界条件
- NEMO 的通用 `BC_INLET`
- 涉及状态变量变化的全链路兼容
- 刚性源项 Jacobian 一致性
- 新增自定义机理后与输运模型、弛豫模型之间的耦合一致性

### 8.3 单温模型的潜在注意事项

尽管单温模型已系统接通，但它本质上是对双温结构的退化实现。因此后续新增任何功能时，都必须明确判断：

- 该功能是仅适用于双温模型，还是单温也应支持。
- 若单温支持，是否应置零某些项，还是应改为合并处理。
- 输出与重启变量是否仍保持一致。

这在添加新源项或新能量模式时尤为重要。

---

## 9. 对当前分支的技术判断

基于当前代码阅读，可以明确判断：

1. 当前仓库中的 NEMO 二次开发是真实存在且规模较大的，不是仅停留在配置层。
2. 单温模型已经贯穿物理、数值与输出主链路。
3. Mutation++ 的局部点隐式化学方法是当前分支的重要定制能力之一。
4. SU2_NONEQ 的 custom chemistry 扩展也是结构性改动，具有继续扩展机理的良好基础。
5. 网格自适应已经通过 MMG 外循环接入当前分支，并针对 NEMO 变量做了适配。

因此，当前分支已经具备继续开展 NEMO 二次开发的基础，但后续应优先在以下原则下推进：

- 保持状态变量定义与热力学接口的一致性。
- 把物理模型修改、源项 Jacobian 修改和时间推进稳定性一起考虑。
- 对边界条件和壁面模型保持审慎，不默认其已完全成熟。
- 在新增机理或新方程时，同步检查输出、重启、自适应和 Jacobian 链路。

---

## 10. 建议的后续工作顺序

建议后续开发按以下顺序组织：

1. 先整理一份 NEMO 控制方程与代码变量的一一对应表。
2. 再整理单温模型与点隐式化学的完整算法流程图。
3. 然后确定下一步二开的目标属于哪一类：
   - 新物理模型
   - 新反应机理
   - 新数值离散
   - 新边界条件
   - 新自适应指标
4. 最后围绕该目标建立最小修改闭环：
   - 配置
   - 变量
   - 物性
   - 源项
   - Jacobian
   - 输出
   - 测试样例

---

## 11. 附录：本次重点阅读文件

### 宏观框架

- `SU2_CFD/src/SU2_CFD.cpp`
- `SU2_CFD/src/solvers/CSolverFactory.cpp`
- `SU2_CFD/src/iteration/CFluidIteration.cpp`

### NEMO 主求解器

- `SU2_CFD/src/solvers/CNEMOEulerSolver.cpp`
- `SU2_CFD/src/solvers/CNEMONSSolver.cpp`
- `SU2_CFD/include/solvers/CNEMOEulerSolver.hpp`

### 变量与状态恢复

- `SU2_CFD/src/variables/CNEMOEulerVariable.cpp`

### 热力学与化学后端

- `SU2_CFD/include/fluid/CNEMOGas.hpp`
- `SU2_CFD/src/fluid/CNEMOGas.cpp`
- `SU2_CFD/src/fluid/CMutationTCLib.cpp`
- `SU2_CFD/src/fluid/CSU2TCLib.cpp`

### NEMO 数值方法

- `SU2_CFD/src/numerics/NEMO/CNEMONumerics.cpp`
- `SU2_CFD/src/numerics/NEMO/NEMO_sources.cpp`
- `SU2_CFD/src/numerics/NEMO/convection/roe.cpp`
- `SU2_CFD/src/numerics/NEMO/convection/ausm_slau.cpp`
- `SU2_CFD/src/numerics/NEMO/convection/msw.cpp`

### 输出与自适应

- `SU2_CFD/src/output/CNEMOCompOutput.cpp`
- `SU2_CFD/src/drivers/CSinglezoneDriver.cpp`
- `SU2_PY/mesh_adaptation.py`
- `SU2_PY/SU2/adap/mmg.py`
- `SU2_PY/SU2/adap/tools.py`

### 配置入口

- `Common/include/option_structure.hpp`
- `Common/src/CConfig.cpp`

