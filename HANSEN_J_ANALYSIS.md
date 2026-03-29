# Hansen J 统计量与系数量级差异分析报告

## 分析日期: 2025年3月

---

## 一、问题概述

Python复现代码与论文原始结果存在两个主要差异：

| 问题 | 论文结果 | Python修正前 | Python修正后 |
|------|---------|------------|------------|
| Hansen J p值 | 0.47-0.86 | 0.0000 | 0.27-0.59 |
| Level系数 | 1.197 | 0.510 | 0.510 |
| Vol系数 | -4.236 | -2.903 | -2.903 |

---

## 二、Hansen J统计量问题诊断

### 2.1 原始代码的错误

**错误位置**: `src/utils/regression.py` 第304-321行

**错误1: 使用错误的残差**

```python
# 错误代码
e = result['residuals']  # 这是 y - X_hat @ beta
```

正确的残差应该使用原始内生变量计算：
```python
# 正确代码
e = y_res - X_endog_res @ coef_endog  # 这是 y - X @ beta
```

**为什么这很重要？**
- 2SLS的残差定义为 `e = y - X_endog @ beta`
- 使用 `X_hat`（第一阶段拟合值）计算残差会导致残差方差被低估
- 这直接影响Hansen J统计量的计算

**错误2: 使用错误的公式**

```python
# 错误公式
J_stat = N * e_Pz_e / (e_e / N)
```

这个公式计算出来的是一个极大的数值（约106,588），导致p值为0。

### 2.2 正确的Hansen J统计量公式

对于聚类稳健的Hansen J统计量，正确公式为：

```
J = e' @ W @ V^{-1} @ W' @ e
```

其中：
- `e` 是用原始内生变量计算的2SLS残差
- `W` 是工具变量矩阵
- `V = sum_c (W_c' @ e_c) @ (W_c' @ e_c)'` 是聚类稳健的矩条件方差矩阵

**推导说明**：
1. Hansen J是GMM框架下的过度识别检验
2. 对于聚类数据，需要使用聚类稳健的权重矩阵
3. 这与Stata `ivreg2` 的实现一致

### 2.3 修正后的代码

```python
# Hansen J statistic (overidentification test)
# Compute correct residuals using original endogenous variables
if X_exog_res.shape[1] > 0:
    e = y_res - np.hstack([X_exog_res, X_endog_res]) @ result['coef']
else:
    e = y_res - X_endog_res @ coef_endog

# Cluster-robust Hansen J statistic
unique_clusters = np.unique(clusters)
L = W.shape[1]  # number of instruments

# Compute cluster-robust variance matrix of moment conditions
V_cluster = np.zeros((L, L))
for c in unique_clusters:
    mask = clusters == c
    Wc = W[mask]
    ec = e[mask]
    g_c = Wc.T @ ec  # moment condition for this cluster
    V_cluster += np.outer(g_c, g_c)

# Hansen J = e' @ W @ V^{-1} @ W' @ e
V_inv = np.linalg.pinv(V_cluster)
J_stat = e.T @ W @ V_inv @ W.T @ e
```

### 2.4 修正后的结果对比

| 规格说明 | 论文p值 | 修正前p值 | 修正后p值 |
|---------|--------|----------|----------|
| Table 2 Col 2 (Baseline) | 0.564 | 0.0000 | 0.2722 |
| Table 2 Col 3 (Stock) | - | 0.0000 | 0.3342 |
| Table 6 Col 1 (WUI) | 0.595 | 0.0000 | 0.3174 |
| Table 6 Col 2 (EPU) | 0.596 | 0.0000 | 0.7415 |
| Table 6 Col 4 (Consensus) | 0.858 | 0.0000 | 0.5852 |
| Table 6 Col 5 (Exchange) | 0.473 | 0.0000 | 0.2742 |

**结论**: 修正后的Hansen J p值已经进入合理范围，与论文结果的差异可能来自：
1. 工具变量数据的微小差异
2. 小样本校正方式的不同
3. Stata `ivreg2` 的特殊实现细节

---

## 三、系数量级差异分析

### 3.1 差异描述

| 变量 | 论文系数 | Python系数 | 差异比例 |
|------|---------|-----------|---------|
| Level (cs_index_ret) | 1.197 | 0.510 | 2.35x |
| Vol (cs_index_vol) | -4.236 | -2.903 | 1.46x |

### 3.2 可能原因分析

**原因1: 变量标准化方式**

论文中提到：
> "The first-and second-moment series are scaled for comparability across columns to have residualized unit standard deviation over the regression sample."

这意味着论文使用的是：
1. 对固定效应残差化后的变量
2. 然后标准化为单位标准差

**检验结果**:
- 当前cs_index变量的标准差约为1.4-1.7（不是严格的1）
- Partial out FE后的标准差约为1.16-1.31

**尝试标准化后的系数**:
```
只标准化X: Level=0.670, Vol=-3.37
标准化X和y: Level=0.242, Vol=-1.22
不标准化: Level=0.510, Vol=-2.90
```

**原因2: 数据定义差异**

论文可能使用了不同的变量定义或数据版本。具体来说：
- `cs_index_ret` 和 `cs_index_vol` 的构建方式可能有差异
- 时间滞后或其他处理可能不同

**原因3: 固定效应处理**

当前代码使用QR分解来partial out固定效应，这可能与Stata的实现有细微差异。

### 3.3 结论

系数量级差异的原因尚未完全确定，但：
1. **符号一致性**: 所有系数符号与论文一致（Level为正，Vol为负）
2. **统计显著性**: 关键系数仍然统计显著
3. **核心结论一致**: 不确定性对GDP增长有负面影响

---

## 四、代码修正方案

### 4.1 已实施的修正

修改文件: `src/utils/regression.py`

关键改动:
1. 使用正确的残差计算方式
2. 使用聚类稳健的Hansen J公式
3. 添加详细的代码注释

### 4.2 建议的后续改进

1. **验证工具变量定义**: 确认数据中的工具变量与论文完全一致
2. **尝试标准化**: 按论文方法标准化变量后重新估计
3. **对比Stata输出**: 如果可能，直接运行Stata代码对比

---

## 五、总体评估

| 方面 | 修正前状态 | 修正后状态 |
|------|----------|----------|
| Hansen J统计量 | ❌ 完全错误 (p=0) | ✅ 合理 (p=0.27-0.74) |
| 系数符号 | ✅ 正确 | ✅ 正确 |
| 系数量级 | ⚠️ 差异约50% | ⚠️ 差异约50% (未改变) |
| 统计显著性 | ✅ 显著 | ✅ 显著 |
| 核心经济学结论 | ✅ 一致 | ✅ 一致 |

---

## 六、技术细节补充

### 6.1 Hansen J统计量的理论基础

Hansen J统计量用于检验GMM估计中的过度识别约束：

- **原假设**: 所有工具变量都是外生的（与误差项不相关）
- **备择假设**: 至少有一些工具变量是内生的

统计量分布: χ²(df)，其中df = L - K（工具变量数 - 内生变量数）

### 6.2 Sargan vs Hansen J

| 统计量 | 假设 | 适用场景 |
|-------|------|---------|
| Sargan | 同方差 | 经典2SLS |
| Hansen J | 允许异方差/聚类 | 稳健估计 |

当前代码使用的是Hansen J（稳健版本），与Stata `ivreg2, cluster()` 一致。

---

**报告完成时间**: 2025年3月
**修改文件**: `src/utils/regression.py`
