# FIRST_STAGE.m <-> Python 对齐清单（逐条）

## 1) DATAMAT 列定义与构造
- MATLAB `FIRST_STAGE.m` 列定义：`22-35`
- Python 对应 DATAMAT 构造：`src/model/gmm.py:166-197`
- 当前状态：已对齐（含 `ct` 与 `ct-1` 滞后写入规则）。

## 2) 一阶矩标准化（总体标准差）
- MATLAB：`DATAMAT(:,4)/std(...,1)` 与 `DATAMAT(:,6)/std(...,1)`（`37-39`）
- Python：`_std_pop(ddof=0)` + 对应缩放（`200-211`）
- 当前状态：已对齐。

## 3) 二阶矩对数化后标准化
- MATLAB：`log` + `/std(...,1)`（`41-46`）
- Python：`log` + `_std_pop(ddof=0)`（`213-221`）
- 当前状态：已对齐。

## 4) dummyvar 处理
- MATLAB：`dummyvar(country)`、`dummyvar(time)`（`48-50`）
- Python：按出现类别生成完整哑变量（`223-233`）
- 当前状态：已对齐（不丢基类）。

## 5) 工具变量矩阵 Z
- MATLAB：`Z=[灾难4列, 国家哑变量, 时间哑变量]`（`52`）
- Python：`Z=np.column_stack([dm[:,7:11], ccv, ttv])`（`235`）
- 当前状态：已对齐。

## 6) 一阶段回归与预测
- MATLAB：宏观一阶段 `55-60`，微观一阶段 `67-73`
- Python：宏观一阶段 `243-248`，微观一阶段 `251-256`
- 当前状态：已对齐。

## 7) 二阶段回归（growthsimyr）
- MATLAB：宏观 `63-65`，微观 `76-78`
- Python：宏观 `247-249`，微观 `255-256`
- 当前状态：已对齐。

## 8) 求解器语义（MATLAB 反斜杠）
- MATLAB：`(X'X)\(X'y)`（多处）
- Python：统一 `np.linalg.lstsq`（`237-240`）
- 当前状态：已对齐到最小二乘语义。

## 9) OUTVEC 映射顺序
- MATLAB：`85-86`
- Python：`258-263`
- 当前状态：已对齐（20 个矩次序一致）。

## 10) uncexogsim 随机过程（Fortran）
- Fortran 抽样顺序：`551-556, 564-573`
- Python 对应：`src/model/gmm.py:302-318` + `src/model/simulation.py` 的外部覆盖输入
- 当前状态：已按顺序消费并接入仿真（含 `ascorrshocks` 规则，`397-400`）。
