# GE Solver & GMM 模块 Fortran vs Python 详细对比

## 1. GE 求解器 - 市场出清

### Fortran (行 1206-1218)
```fortran
! 超额需求
ep0(piter) = 1.0/pval - Cvalp

! 存储加总变量
Cp0(piter) = Cvalp
Yp0(piter) = Yvalp
...

! 市场出清价格和变量
psim(t) = pval
Csim(t) = Cvalp
Ysim(t) = Yvalp
Isim(t) = Ivalp
...
Ksim(t+1) = Kbarprimevalp
```

### Python
```python
# 超额需求
excess_demand = 1.0 / price - C_agg

# 存储结果
p_sim[t] = price
C_sim[t] = C_agg
Y_sim[t] = Y_agg
I_sim[t] = I_agg
...
K_sim[t + 1] = K_prime
```

**差异分析：** ✅ 逻辑一致

## 2. GMM 目标函数

### Fortran (fGMM 函数)
```fortran
double precision function fGMM(x)
implicit none

double precision :: x(8)  ! 8个参数

! 运行完整模型模拟
! 计算20个矩

! GMM 目标
GMMobj = 0.0
do ct=1,nummom
    GMMobj = GMMobj + ((MATLAB_MOMS(ct) - DATA_MOMS(ct))/DATA_SE(ct))**2
end do

fGMM = GMMobj
end function
```

### Python
```python
def gmm_objective(x, params, grids, data_moments, data_se, vfi_sol):
    # 计算20个模拟矩
    sim_moments = compute_simulated_moments(...)
    
    # GMM 目标
    diff = (sim_moments - data_moments) / data_se
    gmm_value = np.sum(diff ** 2)
    
    return gmm_value
```

**差异分析：** ✅ 逻辑一致

## 3. 矩匹配结构

### 20个矩的组织
| 索引 | Fortran | Python | 说明 |
|------|---------|--------|------|
| 1-4 | FIRST_STAGE_MACRO_DATA(:,1) | first_stage_macro[:, 0] | 宏观样本，水平LHS |
| 5-8 | FIRST_STAGE_MACRO_DATA(:,2) | first_stage_macro[:, 1] | 宏观样本，波动LHS |
| 9-10 | SECOND_STAGE_MACRO_DATA | second_stage_macro | 宏观样本，第二阶段 |
| 11-14 | FIRST_STAGE_MICRO_DATA(:,1) | first_stage_micro[:, 0] | 微观样本，水平LHS |
| 15-18 | FIRST_STAGE_MICRO_DATA(:,2) | first_stage_micro[:, 1] | 微观样本，波动LHS |
| 19-20 | SECOND_STAGE_MICRO_DATA | second_stage_micro | 微观样本，第二阶段 |

**差异分析：** ✅ 结构完全一致

## 4. 参数向量

### Fortran 索引 (行 39-47)
```fortran
! 1 - nat_dis levels
! 2 - pol_shock levels
! 3 - revolution levels
! 4 - terrorist levels
! 5 - nat_dis unc
! 6 - pol_shock unc
! 7 - revolution unc
! 8 - terrorist unc
```

### Python
```python
def get_param_vector(self):
    return np.array([
        self.disaster_levels['nat_dis'],      # x[0]
        self.disaster_levels['pol_shock'],    # x[1]
        self.disaster_levels['revolution'],   # x[2]
        self.disaster_levels['terrorist'],    # x[3]
        self.disaster_unc_probs['nat_dis'],   # x[4]
        self.disaster_unc_probs['pol_shock'], # x[5]
        self.disaster_unc_probs['revolution'],# x[6]
        self.disaster_unc_probs['terrorist']  # x[7]
    ])
```

**差异分析：** ✅ 索引对应正确（考虑 0-based vs 1-based）

## 5. 数据矩

### Fortran (行 214-260)
```fortran
FIRST_STAGE_MACRO_DATA(1,:) = (/ -0.071, -0.028 /) !Nat Dis
FIRST_STAGE_MACRO_DATA(2,:) = (/ 1.657, 1.693 /)   !Pol Shock
FIRST_STAGE_MACRO_DATA(3,:) = (/ -6.154, 7.841 /)  !Revolution
FIRST_STAGE_MACRO_DATA(4,:) = (/ -0.047, -0.011 /) !Terrorist

SECOND_STAGE_MACRO_DATA = (/1.557, -3.859/)

FIRST_STAGE_MICRO_DATA(1,:) = (/ -0.147, 0.004/)   !Nat Dis
FIRST_STAGE_MICRO_DATA(2,:) = (/ 1.852, 0.508 /)   !Pol Shock
FIRST_STAGE_MICRO_DATA(3,:) = (/ -4.818, 3.201 /)  !Revolution
FIRST_STAGE_MICRO_DATA(4,:) = (/ -0.117, 0.133 /)  !Terrorist

SECOND_STAGE_MICRO_DATA = (/0.736, -9.735/)
```

### Python (params.py)
```python
first_stage_macro: np.ndarray = field(default_factory=lambda: np.array([
    [-0.071, -0.028],   # Nat Dis
    [1.657, 1.693],     # Pol Shock
    [-6.154, 7.841],    # Revolution
    [-0.047, -0.011]    # Terrorist
]))

second_stage_macro: np.ndarray = field(default_factory=lambda: np.array([1.557, -3.859]))

first_stage_micro: np.ndarray = field(default_factory=lambda: np.array([
    [-0.147, 0.004],    # Nat Dis
    [1.852, 0.508],     # Pol Shock
    [-4.818, 3.201],    # Revolution
    [-0.117, 0.133]     # Terrorist
]))

second_stage_micro: np.ndarray = field(default_factory=lambda: np.array([0.736, -9.735]))
```

**差异分析：** ✅ 数值完全一致

## 6. 标准误

### Fortran (行 269-301)
```fortran
DATA_SE(1) = 0.1059
DATA_SE(2) = 0.0551
DATA_SE(3) = 1.0835
...
DATA_SE(20) = 1.533
```

### Python
```python
data_se: np.ndarray = field(default_factory=lambda: np.array([
    0.1059, 0.0551, 1.0835, 0.0514,  # First stage, levels LHS, macro
    0.0823, 0.1159, 2.2358, 0.0490,  # First stage, vol LHS, macro
    0.2906, 0.2844,                  # Second stage, growth LHS, macro
    0.112, 0.085, 1.198, 0.044,      # First stage, levels LHS, micro
    0.102, 0.130, 1.275, 0.083,      # First stage, vol LHS, micro
    0.558, 1.533                     # Second stage, growth LHS, micro
]))
```

**差异分析：** ✅ 数值完全一致

## 7. PSO 优化参数

### Fortran (行 29-37)
```fortran
npart = 75
xpsotol = 1.0e-3
fpsotol = 1.0e-3
phipso = (/2.05,2.05/)
maxpsoit = 5000
```

### Python (optimizer.py)
```python
n_particles = 75
x_tol = 1e-3
f_tol = 1e-3
phi = (2.05, 2.05)
max_iter = 5000
```

**差异分析：** ✅ 参数完全一致

