# GitHub 项目模块完整性分析报告

## 一、项目概览

### 项目结构
```
uncertainty-disasters/
├── data/
│   ├── IV/
│   │   ├── panel_iv_data.dta (8.9MB) ✅
│   │   └── dstats.csv (1.3KB) ✅
│   ├── IV_VAR/
│   │   └── VARdata.csv (195KB) ✅
│   └── LMN_VAR/
│       └── Dates_and_Data.dta (470KB) ✅
├── src/
│   ├── iv/
│   │   ├── __init__.py (21B)
│   │   ├── __main__.py (131B)
│   │   └── panel_iv.py (17.4KB) ✅
│   ├── iv_var/
│   │   ├── __init__.py (25B)
│   │   ├── __main__.py (127B)
│   │   └── estimation.py (10.8KB) ✅
│   ├── lmn_var/
│   │   ├── __init__.py (26B)
│   │   ├── __main__.py (126B)
│   │   └── estimation.py (6.2KB) ⚠️
│   ├── model/
│   │   ├── __init__.py (24B)
│   │   ├── __main__.py (141B)
│   │   └── solve.py (~3.6KB) ❌
│   └── utils/
│       ├── __init__.py (24B)
│       └── regression.py (8.3KB) ✅
└── README.md
```

---

## 二、模块详细分析

### 2.1 IV 模块 (Tables 1-6)

**代码文件**: `src/iv/panel_iv.py` (17,429 字节)

**论文输出**: Tables 1-6

**数据文件**: `data/IV/panel_iv_data.dta` (8.9MB) ✅

#### 功能完整性检查

| 功能 | 论文要求 | Python 实现 | 状态 |
|------|---------|------------|------|
| **Table 1** | 描述性统计 | `table1()` 方法 | ✅ 完整 |
| **Table 2** | 基准 IV 回归 (5列) | `table2()` 方法 | ✅ 完整 |
| **Table 3** | 稳健性检验 (6列) | `table3()` 方法 | ✅ 完整 |
| **Table 4** | 贸易/距离加权 (6列) | `table4()` 方法 | ✅ 完整 |
| **Table 5** | 媒体权重 (6列) | `table5()` 方法 | ✅ 完整 |
| **Table 6** | 替代不确定性代理 (5列) | `table6()` 方法 | ✅ 完整 |

#### 技术实现检查

| 技术点 | Stata 原始 | Python 复现 | 状态 |
|--------|-----------|------------|------|
| 国家固定效应 | `areg, ab(country)` | `demean_by_group()` | ✅ |
| 时间固定效应 | `i.yq_int` | `demean_by_group()` | ✅ |
| 聚类标准误 | `cluster(country)` | Cluster-robust SE | ✅ |
| IV/2SLS | `ivreg2` | `iv2sls()` | ✅ |
| Hansen J 检验 | `ivreg2` 输出 | 计算 J 统计量 | ✅ |
| 第一阶段 F 统计量 | `ivreg2` 输出 | 计算F统计量 | ✅ |
| 分析权重 | `[aw=lpop]` | 加权 IV | ✅ |
| FWL 定理 | `reghdfe` | 迭代去均值 | ✅ |

#### 潜在问题

1. **数据格式**: 原始使用 `.dta` 文件，Python 使用 `pd.read_stata()` 兼容 ✅
2. **变量名**: 代码中变量名与 Stata 全局变量命名一致 ✅
3. **工具变量集**: 包含完整的 8 组工具变量集 (iv, d_iv, t_iv, s0, sMed, s2, sprd, nws) ✅

**结论**: **IV 模块完整度约 95%**，可正确复现 Tables 1-6。

---

### 2.2 IV_VAR 模块 (Figures 6-7)

**代码文件**: `src/iv_var/estimation.py` (10,809 字节)

**论文输出**: Figures 6-7

**数据文件**: `data/IV_VAR/VARdata.csv` (195KB) ✅

#### 功能完整性检查

| 功能 | 论文要求 | Python 实现 | 状态 |
|------|---------|------------|------|
| VAR 估计 | 3变量 VAR (ydgdp, avgret, lavgvol) | `_build_moments()` | ✅ |
| GMM 目标函数 | 18 个矩条件 | `_gmm_obj()` | ✅ |
| IRF 计算 | 50 期脉冲响应 | `_compute_irf()` | ✅ |
| Bootstrap SE | Stationary block bootstrap | `bootstrap_se()` | ✅ |
| **Figure 6** | 基准 IRF 图 | `plot_fig6()` | ✅ |
| **Figure 7** | 稳健性 IRF 图 | `plot_fig7()` | ✅ |

#### 技术实现检查

| 技术点 | MATLAB 原始 | Python 复现 | 状态 |
|--------|------------|------------|------|
| VAR 滞后 | 3 lags | `Nlags = 3` | ✅ |
| GMM 优化 | `fminsearch` | `minimize(method='Nelder-Mead')` | ✅ |
| 参数数量 | 17 个参数 | `Nparams = 17` | ✅ |
| 矩条件数量 | 18 个矩 | `Nmoms = 18` | ✅ |
| Bootstrap 次数 | 500 | `n_boot=500` | ✅ |
| Block size | 4 | `block_size=4` | ✅ |
| 置信区间 | 90% (1.645) | `1.645*SE` | ✅ |

#### 代码问题

⚠️ **发现 Bug**: `_stationary_bb()` 方法中有语法错误：
```python
# 第 155-165 行存在错误:
Xb, j] = X_ext[I]+jj, j]  # 应为 Xb[h, j] = X_ext[I[m]+jj, j]
Db, j] = D_ext[I]+jj, j]  # 应为 Db[h, j] = D_ext[I[m]+jj, j]
```

这是索引语法错误，会导致 bootstrap 无法运行。

**结论**: **IV_VAR 模块完整度约 85%**，存在代码 Bug 需要修复。

---

### 2.3 LMN_VAR 模块 (Figures 3-5)

**代码文件**: `src/lmn_var/estimation.py` (6,186 字节)

**论文输出**: Figures 3-5

**数据文件**: `data/LMN_VAR/Dates_and_Data.dta` (470KB) ✅

#### 功能完整性检查

| 功能 | 论文要求 | Python 实现 | 状态 |
|------|---------|------------|------|
| VAR + FE 估计 | `reghdfe` | `step1_var_fe()` | ✅ |
| Admissible IRF | 100,000 抽取 | `step2_admissible()` | ✅ |
| 集合识别 | Ludvigson et al. | 中位数 + 68% CI | ✅ |
| **Figure 3** | GDP IRF | `step3_figures()` | ✅ |
| **Figure 4** | Uncertainty IRF | `step3_figures()` | ✅ |
| **Figure 5** | Returns IRF | `step3_figures()` | ✅ |

#### 技术实现检查

| 技术点 | Stata+MATLAB 原始 | Python 复现 | 状态 |
|--------|------------------|------------|------|
| 双固定效应 | country + time | 迭代去均值 | ✅ |
| VAR 滞后 | 3 lags | `Nlags = 3` | ✅ |
| QR 分解 | 正交化冲击 | `np.linalg.qr()` | ✅ |
| 约束条件 | GDP↓, Vol↑ | IRF 约束检查 | ✅ |
| 抽取数量 | 100,000 | `n_draws=100000` | ✅ |
| 分位数 | 16%, 50%, 84% | `np.percentile()` | ✅ |

#### 潜在问题

1. **计算效率**: 100,000 次抽取可能需要优化
2. **内存使用**: 存储 admissible IRFs 可能占用大量内存
3. **简化实现**: 与原始 MATLAB 代码相比有些简化

**结论**: **LMN_VAR 模块完整度约 80%**，基本可用。

---

### 2.4 MODEL 模块 (Figure 8)

**代码文件**: `src/model/solve.py` (~3,600 字节，约 130 行)

**论文输出**: Figure 8

**数据文件**: 无 (模型模拟)

#### 功能完整性检查

| 功能 | 论文要求 | Python 实现 | 状态 |
|------|---------|------------|------|
| 参数校准 | ~50 个参数 | 9 个简化参数 | ❌ 不完整 |
| 状态空间 | 5 维 + 加总 | 2 维 | ❌ 错误 |
| VFI | Howard 加速 + 并行 | 简单迭代 | ❌ 不完整 |
| 调整成本 | 非凸 (资本+劳动) | 二次型 | ❌ 错误 |
| 灾难事件 | 4 种类型 | 无 | ❌ 缺失 |
| GMM 估计 | 20 个矩匹配 | 无 | ❌ 缺失 |
| PSO 优化 | 75 粒子, 5000 迭代 | 无 | ❌ 缺失 |
| IV 回归 | MATLAB 计算 | 无 | ❌ 缺失 |
| **Figure 8** | GDP IRF | 简化版 | ⚠️ 骨架代码 |

#### 详细对比

| 项目 | 原始 Fortran | Python 复现 | 差异 |
|------|-------------|------------|------|
| 代码规模 | ~3,000 行 | ~130 行 | -96% |
| 参数数量 | ~50 | 9 | -82% |
| 状态维度 | 17M 点 | 10K 点 | -99% |
| 企业数量 | 800 | 无 | 缺失 |
| 劳动变量 | 有 | 无 | 缺失 |
| 调整成本类型 | 非凸 | 二次型 | 错误 |

**结论**: **MODEL 模块完整度 < 5%**，需要完全重写。

---

### 2.5 Utils 模块

**代码文件**: `src/utils/regression.py` (8,273 字节)

#### 功能完整性检查

| 函数 | 功能 | 状态 |
|------|------|------|
| `demean_by_group()` | 单固定效应去均值 | ✅ |
| `partial_out_fe()` | 多固定效应去均值 (FWL) | ✅ |
| `ols_with_se()` | OLS + 聚类/稳健 SE | ✅ |
| `iv2sls()` | 两阶段最小二乘 | ✅ |
| `get_cc_yy_cols()` | 提取虚拟变量列 | ✅ |
| `format_coef()` | 格式化输出 | ✅ |

**结论**: **Utils 模块完整度 95%**，功能完整。

---

## 三、数据文件完整性

| 模块 | 所需数据 | 存在 | 大小 | 状态 |
|------|---------|------|------|------|
| IV | panel_iv_data.dta | ✅ | 8.9MB | 完整 |
| IV | dstats.csv | ✅ | 1.3KB | 完整 |
| IV_VAR | VARdata.csv | ✅ | 195KB | 完整 |
| LMN_VAR | Dates_and_Data.dta | ✅ | 470KB | 完整 |
| MODEL | 无需数据 | N/A | N/A | N/A |

**所有数据文件完整** ✅

---

## 四、总体评估

### 完整度汇总

| 模块 | 论文输出 | 代码完整度 | 数据完整度 | 可复现性 |
|------|---------|-----------|-----------|---------|
| **IV** | Tables 1-6 | **95%** ✅ | 100% ✅ | **高** |
| **IV_VAR** | Figures 6-7 | **85%** ⚠️ | 100% ✅ | **中** (有Bug) |
| **LMN_VAR** | Figures 3-5 | **80%** ⚠️ | 100% ✅ | **中** |
| **MODEL** | Figure 8 | **<5%** ❌ | N/A | **低** |

### 关键问题清单

#### 高优先级 (必须修复)

1. **IV_VAR Bootstrap Bug**: `_stationary_bb()` 中的索引语法错误
2. **MODEL 模块缺失**: 需要完全重写，参考原始 Fortran 代码

#### 中优先级 (建议修复)

3. **LMN_VAR 效率**: 可优化大样本抽取的性能
4. **IV_VAR IRF 缩放**: Bootstrap 中的缩放逻辑需验证

#### 低优先级 (可选)

5. **输出格式**: 可添加 LaTeX 表格输出
6. **文档**: 可添加更多代码注释

---

## 五、修正建议

### 5.1 IV_VAR 模块修正

**Bug 位置**: `src/iv_var/estimation.py` 第 155-165 行

**错误代码**:
```python
Xb, j] = X_ext[I]+jj, j]
Db, j] = D_ext[I]+jj, j]
```

**修正代码**:
```python
def _stationary_bb(self, X, D, rng, block_size):
    T, NX = X.shape
    ND = D.shape[1]
    X_ext = np.vstack([X, X[:T-1]])
    D_ext = np.vstack([D, D[:T-1]])
    n_ext = T - 1

    I = np.round(rng.uniform(0, 1, T) * (n_ext - 1)).astype(int)
    b = rng.geometric(1.0 / block_size, T)

    Xb = np.zeros((T, NX))
    Db = np.zeros((T, ND))
    
    for j in range(NX):
        h = 0
        for m in range(T):
            for jj in range(b[m]):
                if h >= T: break
                Xb[h, j] = X_ext[I[m]+jj, j]
                h += 1
            if h >= T: break
    
    for j in range(ND):
        h = 0
        for m in range(T):
            for jj in range(b[m]):
                if h >= T: break
                Db[h, j] = D_ext[I[m]+jj, j]
                h += 1
            if h >= T: break
    
    return Xb, Db
```

### 5.2 MODEL 模块重写

需要按照原始 Fortran 代码重新实现，预计需要：
- 参数模块 (~200 行)
- 网格生成 (~300 行)
- 数值工具 (~400 行)
- 调整成本 (~200 行)
- 价值函数迭代 (~500 行)
- 模拟系统 (~700 行)
- 估计系统 (~500 行)
- PSO 优化 (~300 行)

**预计总代码量**: ~3,300 行

---

## 六、结论

### 可以复现的论文内容

| 内容 | 状态 |
|------|------|
| Table 1 (描述性统计) | ✅ 可复现 |
| Table 2 (基准结果) | ✅ 可复现 |
| Table 3 (稳健性检验) | ✅ 可复现 |
| Table 4 (贸易/距离加权) | ✅ 可复现 |
| Table 5 (媒体权重) | ✅ 可复现 |
| Table 6 (替代不确定性代理) | ✅ 可复现 |
| Figure 3 (GDP IRF) | ⚠️ 基本可复现 |
| Figure 4 (Uncertainty IRF) | ⚠️ 基本可复现 |
| Figure 5 (Returns IRF) | ⚠️ 基本可复现 |
| Figure 6 (基准 IRF) | ⚠️ 需修复 Bug |
| Figure 7 (稳健性 IRF) | ⚠️ 需修复 Bug |
| Figure 8 (模型 IRF) | ❌ 无法复现 |

### 总体完成度

- **IV 模块**: 95% ✅ (Tables 1-6 可完整复现)
- **IV_VAR 模块**: 85% ⚠️ (存在 Bug 需修复)
- **LMN_VAR 模块**: 80% ⚠️ (基本可用，可优化)
- **MODEL 模块**: <5% ❌ (需要完全重写)

**项目整体完成度**: 约 **70%**

**建议**: 
1. 优先修复 IV_VAR 的 Bootstrap Bug
2. 重点重写 MODEL 模块 (工作量最大)
3. 优化 LMN_VAR 的计算效率
