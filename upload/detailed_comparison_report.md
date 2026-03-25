# 代码功能逐一检验对比报告

## 一、检验概要

**检验日期**: 2025年9月
**检验方法**: 逐功能对比 Python 复现代码与原始 Stata/MATLAB/Fortran 代码
**目的**: 找出还需完善的部分

---

## 二、IV 模块详细对比

### 2.1 Table 1: 描述性统计

| 项目 | 原始 Stata | Python 复现 | 状态 |
|------|-----------|------------|------|
| 样本筛选 | `if ydgdp~=.` | `self.df['ydgdp'].notna()` | ✅ 正确 |
| 统计量 | count, mean, p50, sd, min, max | 相同 | ✅ 正确 |
| 输出格式 | CSV | CSV | ✅ 正确 |

**结论**: Table 1 完整复现 ✅

---

### 2.2 Table 2: 基准回归

| 项目 | 原始 Stata | Python 复现 | 状态 | 问题 |
|------|-----------|------------|------|------|
| Col 1 OLS | `areg ... ab(country)` | `_run_areg()` | ✅ | |
| Col 2 IV | `ivreg2 ... partial(yy* cc*)` | `_run_iv()` | ✅ | |
| Col 3 IV | `l1avgret, l1lavgvol` | 相同 | ✅ | |
| Col 4 IV | common sample | `sample_filter` | ✅ | |
| Col 5 IV | cross-section | 相同 | ✅ | |
| Hansen J | `e(jp)` | `J_pval` | ✅ | |
| F-stat | `e(F)` | `first_stage[].F_stat` | ✅ | |

**潜在问题**:

1. **聚类标准误计算差异**
   - 原始 Stata: `cluster(country)` 使用 `ivreg2` 的实现
   - Python: 自定义 `ols_with_cluster_se()` 和 `iv2sls_with_cluster_se()`
   - **需要验证**: 自由度校正公式是否完全一致

2. **固定效应偏出方法**
   - 原始 Stata: `partial(yy* cc*)` 使用 `ivreg2` 内部实现
   - Python: QR 分解方法
   ```python
   # src/utils/regression.py 第 212-228 行
   Q, R = np.linalg.qr(P, mode='reduced')
   rank = np.sum(diag_R > 1e-10 * diag_R[0])
   Q_rank = Q[:, :rank]
   Px = Q_rank @ Q_rank.T
   ```
   - **需要验证**: 与 Stata 的数值精度差异

---

### 2.3 Table 3: 稳健性检验

| 项目 | 原始 Stata | Python 复现 | 状态 |
|------|-----------|------------|------|
| Col 1 基准 | ✅ | ✅ | ✅ |
| Col 2 人口加权 | `[aw=lpop]` | `_run_iv_weighted()` | ⚠️ |
| Col 3 加偏度 | ✅ | ✅ | ✅ |
| Col 4 仅波动率和偏度 | ✅ | ✅ | ✅ |
| Col 5 HAR 调整(日度) | ✅ | ✅ | ✅ |
| Col 6 HAR 调整(季度) | ✅ | ✅ | ✅ |

**问题 - 人口加权 IV**:

原始 Stata:
```stata
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol= $iv) [aw=lpop], cluster(country) partial(yy* cc*)
```

Python 实现:
```python
# 第 507 行
w = np.sqrt(data[weight_var].values.astype(np.float64))  # analytic weights
y_w = y * w
X_endog_w = X_endog * w[:, np.newaxis]
```

**潜在问题**: 
- Stata 的 `aw` (analytical weights) 是在回归前对数据进行加权
- Python 实现看似正确，但需要验证权重处理是否完全匹配 Stata

---

### 2.4 Table 4: 贸易/距离加权工具变量

| 项目 | 原始 Stata | Python 复现 | 状态 |
|------|-----------|------------|------|
| 贸易加权工具变量 | `$t_iv` | `self.t_iv` | ✅ |
| 距离加权工具变量 | `$d_iv` | `self.d_iv` | ✅ |
| 公共样本 | `if sample1==1` | `sample_filter=common_mask` | ✅ |

**结论**: Table 4 完整复现 ✅

---

### 2.5 Table 5: 媒体权重

| 项目 | 原始 Stata | Python 复现 | 状态 |
|------|-----------|------------|------|
| 基准 | `$iv` | `self.iv` | ✅ |
| 未缩放 | `l1s0avg...` | ✅ | ✅ |
| 超中位数跳跃 | `l1sMedavg...` | ✅ | ✅ |
| 窄窗口 | `l1s2avg...` | ✅ | ✅ |
| 缩放跳跃回归 | `l1sprdavg...` | ✅ | ✅ |
| 非西方 | `l1nwsavg...` | ✅ | ✅ |

**结论**: Table 5 完整复现 ✅

---

### 2.6 Table 6: 替代不确定性代理

| 项目 | 原始 Stata | Python 复现 | 状态 | 问题 |
|------|-----------|------------|------|------|
| WUI | `l1lavgWUI` | ✅ | ✅ | |
| EPU | `l1lavgEPU` | ✅ | ✅ | |
| EPU+WUI | 合并变量 | `epu_wui` | ⚠️ | 见下 |
| 共识预测 | `l1lgdp_for_sd` | ✅ | ✅ | |
| 汇率 | `l1lavgexchgvol` | ✅ | ✅ | |

**问题 - EPU+WUI 合并变量**:

原始 Stata:
```stata
egen any_epu = mean(EPU), by(country)
gen epu_wui = l1lavgEPU
replace epu_wui = l1lavgWUI if l1lavgEPU==. & any_epu==.
```

Python:
```python
# 第 682-686 行
any_epu = data.groupby('country')['EPU'].transform('mean')
epu_wui = data['l1lavgEPU'].copy()
wui_mask = data['l1lavgEPU'].isna() & any_epu.isna()
epu_wui[wui_mask] = data.loc[wui_mask, 'l1lavgWUI']
```

**问题**: 
- 原始条件是 `l1lavgEPU==. & any_epu==.`，但 `any_epu` 是 EPU 的组均值，不是缺失值
- Python 条件 `any_epu.isna()` 可能与原始意图不符
- **建议修正**: 检查 `any_epu` 是否为缺失值（表示该国家完全没有 EPU 数据）

---

### 2.7 IV 模块总结

| 功能 | 状态 | 备注 |
|------|------|------|
| Table 1 | ✅ 完整 | |
| Table 2 | ✅ 完整 | 需验证聚类 SE 公式 |
| Table 3 | ⚠️ 基本完整 | 人口加权需验证 |
| Table 4 | ✅ 完整 | |
| Table 5 | ✅ 完整 | |
| Table 6 | ⚠️ 基本完整 | EPU+WUI 合并逻辑需修正 |

---

## 三、IV_VAR 模块详细对比

### 3.1 VAR 估计

| 项目 | 原始 MATLAB | Python 复现 | 状态 |
|------|------------|------------|------|
| 变量数 | NX=3 | NX=3 | ✅ |
| 滞后数 | Nlags=3 | Nlags=3 | ✅ |
| 工具变量数 | ND=4 | ND=4 | ✅ |
| 参数数 | Nparams=17 | Nparams=17 | ✅ |
| 矩条件数 | Nmoms=18 | Nmoms=18 | ✅ |

---

### 3.2 GMM 目标函数

原始 MATLAB (fGMMobj.m):
```matlab
% 提取 B 矩阵 (3x3)
B = zeros(NX,NX);
B(1,1) = x(1); B(2,1) = x(2); B(3,1) = x(3);
...

% 计算隐含 Lambda 矩阵
LAMBDA = [1, 0, 0;
          0, sum(Dcoeff(:,1).^2)+1, sum(Dcoeff(:,1).*Dcoeff(:,2));
          0, sum(Dcoeff(:,1).*Dcoeff(:,2)), sum(Dcoeff(:,2).^2)+1];
```

Python:
```python
# 第 140-158 行
B = np.zeros((NX, NX))
B[0, 0] = x[0]; B[1, 0] = x[1]; B[2, 0] = x[2]
...
LAMBDA = np.array([
    [1, 0, 0],
    [0, np.sum(Dcoeff[:, 0] ** 2) + 1, np.sum(Dcoeff[:, 0] * Dcoeff[:, 1])],
    [0, np.sum(Dcoeff[:, 0] * Dcoeff[:, 1]), np.sum(Dcoeff[:, 1] ** 2) + 1],
])
```

**结论**: GMM 目标函数完全匹配 ✅

---

### 3.3 优化器

| 项目 | 原始 MATLAB | Python 复现 | 状态 |
|------|------------|------------|------|
| 方法 | `fminsearch` | `minimize(L-BFGS-B)` | ⚠️ 不同 |
| 最大迭代 | 50000 | 50000 | ✅ |
| 容差 | 默认 | ftol=1e-15, gtol=1e-10 | ⚠️ 可能不同 |

**潜在问题**:
- MATLAB `fminsearch` 使用 Nelder-Mead 单纯形法
- Python 使用 L-BFGS-B 拟牛顿法
- 两种方法可能收敛到略微不同的解

**建议**: 改为使用 Nelder-Mead 以保持一致：
```python
result = minimize(
    self._gmm_objective,
    param0,
    args=(MOMvec,),
    method='Nelder-Mead',  # 匹配 MATLAB fminsearch
    options={'maxiter': 50000, 'xatol': 1e-15, 'fatol': 1e-15},
)
```

---

### 3.4 Bootstrap

| 项目 | 原始 MATLAB | Python 复现 | 状态 |
|------|------------|------------|------|
| 方法 | Stationary block bootstrap | `_stationary_block_bootstrap()` | ✅ |
| 重复次数 | 500 | 500 | ✅ |
| 块大小 | 4 | 4 | ✅ |
| 随机种子 | 3991 | 3991 | ✅ |
| 几何分布 | `sim=1` | `rng.geometric()` | ✅ |

**已修复的 Bug**: 数组索引语法错误已修正 ✅

---

### 3.5 IRF 计算

| 项目 | 原始 MATLAB | Python 复现 | 状态 |
|------|------------|------------|------|
| Companion form | 正确构建 | 正确构建 | ✅ |
| IRF 期数 | lengthIRF=50 | lengthIRF=50 | ✅ |
| 缩放因子 | `sqrt(var(X))*IRF/Bhat` | 相同 | ⚠️ |

**潜在问题 - Bootstrap 缩放**:
```python
# 第 428-434 行
SCALEFACT = (np.sqrt(np.var(X[:, 2])) /
             np.sqrt(np.var(Xb[:, 2])) *
             baseline_result['Bhat'][2, 2] / Bhat_b[2, 2])
```

**问题**: 缺少安全检查，当 `Bhat_b[2,2]=0` 或 `var(Xb[:,2])=0` 时会出错

**建议修正**:
```python
if Bhat_b[2, 2] != 0 and np.var(Xb[:, 2]) > 0:
    SCALEFACT = ...
else:
    boot_bad += 1
    continue
```

---

### 3.6 IV_VAR 模块总结

| 功能 | 状态 | 问题 |
|------|------|------|
| VAR 估计 | ✅ | |
| GMM 目标函数 | ✅ | |
| 优化器 | ⚠️ | 方法不同，建议改用 Nelder-Mead |
| Bootstrap | ✅ | Bug 已修复 |
| 缩放因子 | ⚠️ | 缺少安全检查 |
| IRF 计算 | ✅ | |
| Figure 6-7 | ✅ | |

---

## 四、LMN_VAR 模块详细对比

### 4.1 VAR 固定效应估计

原始 Stata:
```stata
reghdfe depvar l1depvar l2depvar l3depvar l1indepvar1 l1indepvar2 ..., absorb(country yq_int)
```

Python:
```python
# 第 98-106 行
for fe_col in ['country', 'yq_int']:
    groups = eq_data[fe_col].values
    unique_groups = np.unique(groups)
    for g in unique_groups:
        mask = groups == g
        y[mask] -= y[mask].mean()
        X[mask] -= X[mask].mean(axis=0)
```

**结论**: FWL 定理实现正确 ✅

---

### 4.2 变量名映射

| 项目 | 原始 Stata | Python 复现 | 状态 |
|------|-----------|------------|------|
| GDP 增长 | ydgdp | ydgdp | ✅ |
| 一阶矩 | avgret | ret | ⚠️ 不同 |
| 二阶矩 | lavgvol | vol | ⚠️ 不同 |

**问题**: 变量名不一致，需要检查数据文件中的实际列名

---

### 4.3 Admissible Sets 计算

原始 MATLAB (STEP2_MATLAB_ESTIMATION.m):
```matlab
% QR 分解生成正交矩阵
Z = randn(NX,NX);
[Q,R] = qr(Z);
Q = Q * diag(sign(diag(R)));  % 确保 det(Q) = +1

% 结构矩阵
B = chol(Sigma) * Q;
```

Python:
```python
# 第 221-227 行
Z = rng.randn(NX, NX)
Q, R = np.linalg.qr(Z)
Q = Q @ np.diag(np.sign(np.diag(R)))  # 确保 det(Q) = +1

B = Sigma_chol @ Q
```

**结论**: QR 分解实现正确 ✅

---

### 4.4 约束条件

原始 MATLAB 约束（推测）:
- GDP 对不确定性冲击响应为负
- 波动率对不确定性冲击响应为正
- 收益率对不确定性冲击响应为负
- 可能还有其他约束（如持续时间、幅度）

Python 实现:
```python
# 第 271-276 行
restrictions = {
    'gdp_response_at_shock': 'negative',
    'vol_response_at_shock': 'positive',
    'ret_response_at_shock': 'negative',
}
```

**问题**: 约束条件可能过于简化

**建议**: 需要查看原始 MATLAB 代码确认完整约束条件

---

### 4.5 LMN_VAR 模块总结

| 功能 | 状态 | 问题 |
|------|------|------|
| VAR FE 估计 | ✅ | |
| 变量名映射 | ⚠️ | ret/vol vs avgret/lavgvol |
| QR 分解 | ✅ | |
| Cholesky 分解 | ✅ | |
| 约束条件 | ⚠️ | 可能过于简化 |
| IRF 计算 | ✅ | |
| Figure 3-5 | ✅ | |

---

## 五、MODEL 模块详细对比

### 5.1 参数对比

| 参数 | 原始 Fortran | Python 复现 | 状态 |
|------|-------------|------------|------|
| 资本份额 α | 0.25 | 0.33 | ❌ 错误 |
| 劳动份额 ν | 0.5 | 无 | ❌ 缺失 |
| 调整成本参数 θ | 2.0 | 0.5 | ❌ 错误 |
| 资本折旧 δk | 0.026 | 0.025 | ⚠️ 接近 |
| 劳动折旧 δn | 0.088 | 无 | ❌ 缺失 |
| 贴现因子 β | 0.95^0.25 | 0.99 | ❌ 错误 |
| 资本不可逆性 capirrev | 0.339 | 无 | ❌ 缺失 |
| 雇佣成本 hirelin | 0.018*4 | 无 | ❌ 缺失 |
| 解雇成本 firelin | 0.018*4 | 无 | ❌ 缺失 |
| 固定劳动成本 labfix | 0.024*4 | 无 | ❌ 缺失 |

---

### 5.2 网格规模对比

| 项目 | 原始 Fortran | Python 复现 | 状态 |
|------|-------------|------------|------|
| 资本网格 | 150 点 | 500 点 | ⚠️ 不同 |
| 劳动网格 | 75 点 | 无 | ❌ 缺失 |
| 个体生产率网格 | 9 点 | 无 | ❌ 缺失 |
| 总生产率网格 | 21 点 | 无 | ❌ 缺失 |
| 波动率网格 | 2 点 | 20 点 | ⚠️ 不同 |
| 企业数 | 800 | 无 | ❌ 缺失 |

---

### 5.3 调整成本函数对比

原始 Fortran (非凸调整成本):
```fortran
! 资本调整成本
ACk = capirrev * abs(k_new - k_old) + capfix * indicator(k_new != k_old)

! 劳动调整成本
ACl = hirelin * max(l_new - l_old, 0) + firelin * max(l_old - l_new, 0) + labfix * indicator(l_new != l_old)
```

Python (二次型调整成本):
```python
# 第 82-89 行
def _adjustment_cost(self, k_new, k_old):
    p = self.params
    return p['phi'] * (k_new - k_old) ** 2 / (2 * k_old)  # 这是二次型！
```

**严重问题**: 
- 原始模型使用**非凸调整成本**（固定成本 + 可变成本）
- Python 使用**二次型调整成本**
- 这是论文核心机制的错误实现！

---

### 5.4 价值函数迭代对比

原始 Fortran:
```fortran
! 使用 Howard 加速 + 并行化
!$omp parallel private(...)
do accelct=1,accelmaxit
    do ct=1,numstates
        ! 计算 V_new
    end do
end do
```

Python:
```python
# 简单迭代，无 Howard 加速
for iteration in range(max_iter):
    for i in range(Nsigma):
        for j in range(Nk):
            rhs = self._bellman_rhs(V, k, sigma, K_agg)
            best_idx = np.argmax(rhs)
```

**问题**:
1. 无 Howard 加速，收敛速度慢
2. 无并行化
3. 无劳动决策变量
4. 状态空间维度严重不足

---

### 5.5 缺失功能

| 功能 | 原始 Fortran | Python 复现 | 状态 |
|------|-------------|------------|------|
| 灾难事件模拟 | 4 种类型 | 无 | ❌ 缺失 |
| GMM 估计 | 20 个矩匹配 | 无 | ❌ 缺失 |
| PSO 优化 | 75 粒子, 5000 迭代 | 无 | ❌ 缺失 |
| 企业模拟 | 800 家企业 | 无 | ❌ 缺失 |
| MATLAB 接口 | 调用 MATLAB | 无 | ❌ 缺失 |
| 一般均衡 | 价格出清 | 无 | ❌ 缺失 |

---

### 5.6 MODEL 模块总结

| 功能 | 状态 | 问题 |
|------|------|------|
| 参数设置 | ❌ | 大量参数缺失或错误 |
| 状态空间 | ❌ | 维度严重不足 |
| 调整成本 | ❌ | 核心错误：二次型 vs 非凸 |
| 劳动决策 | ❌ | 完全缺失 |
| 价值函数迭代 | ⚠️ | 无 Howard 加速 |
| 灾难事件 | ❌ | 缺失 |
| GMM 估计 | ❌ | 缺失 |
| PSO 优化 | ❌ | 缺失 |
| 企业模拟 | ❌ | 缺失 |
| 一般均衡 | ❌ | 缺失 |

**MODEL 模块完整度: < 5%** ❌

---

## 六、完善建议汇总

### 6.1 高优先级（必须修复）

#### IV 模块
1. **EPU+WUI 合并逻辑修正**
   ```python
   # 当前 (错误):
   wui_mask = data['l1lavgEPU'].isna() & any_epu.isna()
   
   # 修正为:
   wui_mask = data['l1lavgEPU'].isna() & (any_epu == 0)
   ```

#### IV_VAR 模块
2. **优化器方法改为 Nelder-Mead**
   ```python
   result = minimize(
       self._gmm_objective,
       param0,
       args=(MOMvec,),
       method='Nelder-Mead',
       options={'maxiter': 50000, 'xatol': 1e-15, 'fatol': 1e-15},
   )
   ```

3. **Bootstrap 缩放安全检查**
   ```python
   if Bhat_b[2, 2] != 0 and np.var(Xb[:, 2]) > 0:
       SCALEFACT = (np.sqrt(np.var(X[:, 2])) /
                    np.sqrt(np.var(Xb[:, 2])) *
                    baseline_result['Bhat'][2, 2] / Bhat_b[2, 2])
   else:
       boot_bad += 1
       continue
   ```

#### LMN_VAR 模块
4. **变量名验证**
   - 检查 `data/LMN_VAR/Dates_and_Data.dta` 中变量名
   - 如有需要，修改为 `avgret` 和 `lavgvol`

5. **约束条件完善**
   - 查看原始 MATLAB 代码 `STEP2_MATLAB_ESTIMATION.m`
   - 添加完整的约束条件

#### MODEL 模块
6. **完全重写** - 预计需要 ~3,300 行代码

---

### 6.2 中优先级（建议修复）

#### IV 模块
7. **验证聚类标准误公式** - 与 Stata `ivreg2` 对比

#### IV_VAR 模块
8. **验证人口加权 IV** - 与 Stata `[aw=lpop]` 对比

---

### 6.3 低优先级（可选）

9. **添加 LaTeX 表格输出**
10. **添加更多代码注释**
11. **优化计算效率**

---

## 七、结论

### 各模块完善度

| 模块 | 当前状态 | 修复后预期 |
|------|---------|-----------|
| IV | 95% | 98% |
| IV_VAR | 90% | 95% |
| LMN_VAR | 85% | 92% |
| MODEL | < 5% | 需完全重写 |

### 关键问题清单

1. ❗ **MODEL 调整成本函数错误** - 论文核心机制
2. ⚠️ IV_VAR 优化器方法不同
3. ⚠️ LMN_VAR 变量名和约束条件需验证
4. ⚠️ IV EPU+WUI 合并逻辑有误

### 下一步行动

1. 立即修复 IV 和 IV_VAR 模块的小问题
2. 验证 LMN_VAR 变量名和约束条件
3. 着手重写 MODEL 模块（工作量最大）
