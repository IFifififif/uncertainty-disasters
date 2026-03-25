# 项目依赖模块完整清单

## 一、项目信息

**项目名称**: uncertainty-disasters
**项目路径**: `/home/z/my-project/uncertainty-disasters`
**论文**: "Using Disasters to Estimate the Impact of Uncertainty"

---

## 二、直接依赖 (requirements.txt)

### 2.1 核心依赖

| 包名 | 版本要求 | 用途 |
|------|---------|------|
| **pandas** | >=2.0 | 数据处理、读取 Stata/CSV 文件 |
| **numpy** | >=1.24 | 数值计算、矩阵运算 |
| **scipy** | >=1.10 | 优化、插值、统计分布 |
| **statsmodels** | >=0.14 | 统计建模（未直接使用，预留） |
| **linearmodels** | >=5.0 | 面板数据回归（未直接使用，预留） |
| **matplotlib** | >=3.7 | 图表生成 |
| **openpyxl** | >=3.1 | Excel 文件读写 |

---

## 三、模块级依赖分析

### 3.1 IV 模块 (src/iv/)

**功能**: Tables 1-6 面板 IV 回归

**直接依赖**:
```
pandas (pd.read_stata, DataFrame操作)
numpy (np.linalg, 数组运算)
sys, os, pathlib (文件路径)
```

**内部依赖**:
```
src.utils.regression (ols_with_cluster_se, iv2sls_with_cluster_se, demean_multiple_fe)
```

**依赖关系图**:
```
IV 模块
├── pandas (数据处理)
├── numpy (数值计算)
└── utils.regression (回归工具)
    ├── numpy (矩阵运算)
    └── scipy.stats.chi2 (统计检验)
```

---

### 3.2 IV_VAR 模块 (src/iv_var/)

**功能**: Figures 6-7 IV-VAR 估计

**直接依赖**:
```python
import sys
import numpy as np                           # 数值计算
import pandas as pd                          # 数据读取
from pathlib import Path                     # 路径处理
from scipy.optimize import minimize          # GMM 优化
from typing import Optional                  # 类型注解
import matplotlib                            # 图表生成
import matplotlib.pyplot as plt
```

**依赖关系图**:
```
IV_VAR 模块
├── pandas (读取 VARdata.csv)
├── numpy (矩阵运算、随机数)
├── scipy.optimize.minimize (GMM 优化 - L-BFGS-B)
├── scipy.stats (隐式，通过 regression 模块)
└── matplotlib (Figure 6-7 生成)
```

---

### 3.3 LMN_VAR 模块 (src/lmn_var/)

**功能**: Figures 3-5 集合识别 VAR

**直接依赖**:
```python
import sys
import numpy as np                           # 数值计算、QR分解
import pandas as pd                          # 读取 Stata 文件
from pathlib import Path                     # 路径处理
from typing import Optional                  # 类型注解
import matplotlib                            # 图表生成
import matplotlib.pyplot as plt
```

**依赖关系图**:
```
LMN_VAR 模块
├── pandas (读取 Dates_and_Data.dta)
├── numpy (QR分解、矩阵运算、百分位数)
└── matplotlib (Figure 3-5 生成)
```

---

### 3.4 MODEL 模块 (src/model/)

**功能**: Figure 8 结构模型模拟

**直接依赖**:
```python
import sys
import numpy as np                           # 数值计算
from pathlib import Path                     # 路径处理
from scipy.optimize import minimize          # 优化（未完整实现）
from scipy.interpolate import interp1d       # 插值
import matplotlib                            # 图表生成
import matplotlib.pyplot as plt
```

**依赖关系图**:
```
MODEL 模块
├── numpy (矩阵运算、随机数)
├── scipy.optimize (优化 - 部分使用)
├── scipy.interpolate (插值)
└── matplotlib (Figure 8 生成)
```

---

### 3.5 Utils 模块 (src/utils/)

**功能**: 共享回归工具

**直接依赖**:
```python
import numpy as np
import pandas as pd
from typing import Optional
from scipy.stats import chi2                 # Hansen J 检验
```

**提供的函数**:
```
- demean_by_group()           # 单固定效应去均值
- demean_multiple_fe()        # 多固定效应去均值 (FWL)
- ols_with_cluster_se()       # OLS + 聚类标准误
- iv2sls_with_cluster_se()    # 2SLS + 聚类标准误 + Hansen J
- get_cc_yy_cols()            # 提取虚拟变量列
- format_coef_table()         # 格式化输出表格
```

---

## 四、Python 标准库依赖

| 模块 | 用途 | 使用位置 |
|------|------|---------|
| `sys` | 路径处理、命令行参数 | 所有模块 |
| `os` | 文件操作 | IV 模块 |
| `pathlib.Path` | 路径处理 | 所有模块 |
| `typing.Optional` | 类型注解 | 所有模块 |

---

## 五、依赖前置关系图

### 5.1 完整依赖树

```
uncertainty-disasters
│
├── IV 模块 (Tables 1-6)
│   ├── pandas >= 2.0
│   │   └── numpy >= 1.24
│   ├── numpy >= 1.24
│   │   └── (无前置)
│   └── utils.regression
│       ├── numpy
│       └── scipy.stats
│           └── scipy >= 1.10
│
├── IV_VAR 模块 (Figures 6-7)
│   ├── pandas >= 2.0
│   │   └── numpy >= 1.24
│   ├── numpy >= 1.24
│   ├── scipy.optimize.minimize
│   │   └── scipy >= 1.10
│   └── matplotlib >= 3.7
│       ├── numpy
│       └── pillow (可选，图像保存)
│
├── LMN_VAR 模块 (Figures 3-5)
│   ├── pandas >= 2.0
│   │   └── numpy >= 1.24
│   │   └── pyreadstat (读取 .dta 文件)
│   ├── numpy >= 1.24
│   └── matplotlib >= 3.7
│
├── MODEL 模块 (Figure 8)
│   ├── numpy >= 1.24
│   ├── scipy.optimize
│   ├── scipy.interpolate
│   └── matplotlib >= 3.7
│
└── Utils 模块
    ├── numpy >= 1.24
    ├── pandas >= 2.0
    └── scipy.stats
```

---

### 5.2 安装顺序

```bash
# 1. 核心数值计算 (无前置依赖)
pip install numpy>=1.24

# 2. 科学计算 (依赖 numpy)
pip install scipy>=1.10

# 3. 数据处理 (依赖 numpy)
pip install pandas>=2.0

# 4. 统计建模 (依赖 numpy, scipy, pandas)
pip install statsmodels>=0.14
pip install linearmodels>=5.0

# 5. 可视化 (依赖 numpy)
pip install matplotlib>=3.7

# 6. 文件格式支持
pip install openpyxl>=3.1
```

---

## 六、各模块功能-依赖映射

### 6.1 IV 模块功能依赖

| 功能 | 依赖包 | 具体用法 |
|------|--------|---------|
| 读取 Stata 数据 | pandas | `pd.read_stata()` |
| 固定效应去均值 | numpy | 矩阵运算 |
| OLS 回归 | numpy | `np.linalg.solve()` |
| IV/2SLS | numpy | 两阶段最小二乘 |
| 聚类标准误 | numpy | 矩阵运算 |
| Hansen J 检验 | scipy.stats | `chi2.cdf()` |

### 6.2 IV_VAR 模块功能依赖

| 功能 | 依赖包 | 具体用法 |
|------|--------|---------|
| 读取 CSV 数据 | pandas | `pd.read_csv()` |
| VAR 估计 | numpy | 矩阵运算 |
| GMM 优化 | scipy.optimize | `minimize(method='L-BFGS-B')` |
| Bootstrap | numpy | `np.random.RandomState`, `geometric()` |
| IRF 计算 | numpy | `np.linalg.matrix_power()` |
| 图表生成 | matplotlib | `plt.plot()`, `fig.savefig()` |

### 6.3 LMN_VAR 模块功能依赖

| 功能 | 依赖包 | 具体用法 |
|------|--------|---------|
| 读取 Stata 数据 | pandas | `pd.read_stata()` |
| VAR + FE 估计 | numpy | 矩阵运算 |
| QR 分解 | numpy | `np.linalg.qr()` |
| Cholesky 分解 | numpy | `np.linalg.cholesky()` |
| Admissible Sets | numpy | 随机数、矩阵运算 |
| 分位数计算 | numpy | `np.percentile()` |
| 图表生成 | matplotlib | `plt.plot()`, `fill_between()` |

### 6.4 MODEL 模块功能依赖

| 功能 | 依赖包 | 具体用法 |
|------|--------|---------|
| 网格生成 | numpy | `np.linspace()`, `np.exp()` |
| 价值函数迭代 | numpy | 矩阵运算 |
| 插值 | scipy.interpolate | `interp1d()` |
| 优化 | scipy.optimize | `minimize()` |
| 图表生成 | matplotlib | `plt.subplots()` |

---

## 七、依赖版本兼容性

### 7.1 已知兼容版本

| 包 | 最低版本 | 推荐版本 | 备注 |
|---|---------|---------|------|
| Python | 3.9 | 3.11+ | 需要 dataclasses 支持 |
| numpy | 1.24 | 1.26+ | 类型注解改进 |
| pandas | 2.0 | 2.1+ | Stata 读取优化 |
| scipy | 1.10 | 1.12+ | L-BFGS-B 优化改进 |
| matplotlib | 3.7 | 3.8+ | PDF 输出优化 |

### 7.2 潜在兼容性问题

1. **pandas 2.0+**: 读取 Stata 文件需要 `pyreadstat` 包
   ```bash
   pip install pyreadstat
   ```

2. **linearmodels**: 预留但未直接使用，如有需要需安装完整依赖
   ```bash
   pip install linearmodels>=5.0
   # 依赖: statsmodels, scipy, pandas, numpy
   ```

---

## 八、完整安装命令

### 8.1 最小安装 (仅运行项目)

```bash
pip install numpy>=1.24 pandas>=2.0 scipy>=1.10 matplotlib>=3.7 openpyxl>=3.1 pyreadstat
```

### 8.2 完整安装 (包含所有声明依赖)

```bash
pip install -r requirements.txt
pip install pyreadstat  # pandas 读取 Stata 需要
```

### 8.3 开发安装 (包含测试)

```bash
pip install -r requirements.txt
pip install pytest pyreadstat
```

---

## 九、依赖大小估算

| 包 | 大小 (MB) | 安装后 (MB) |
|---|----------|------------|
| numpy | ~25 | ~60 |
| pandas | ~15 | ~100 |
| scipy | ~40 | ~120 |
| matplotlib | ~30 | ~80 |
| statsmodels | ~15 | ~50 |
| linearmodels | ~5 | ~20 |
| openpyxl | ~1 | ~5 |
| pyreadstat | ~2 | ~10 |
| **总计** | ~133 | ~445 |

---

## 十、总结

### 项目依赖层级

```
Level 0 (无前置): numpy
Level 1 (依赖 Level 0): scipy, pandas
Level 2 (依赖 Level 1): statsmodels, matplotlib, linearmodels
Level 3 (依赖 Level 2): 项目模块
```

### 关键依赖

1. **numpy**: 所有数值计算的基础
2. **scipy**: GMM 优化、统计检验
3. **pandas**: 数据读取 (Stata/CSV)
4. **matplotlib**: 所有图表生成

### 可选依赖

1. **statsmodels/linearmodels**: 已声明但未直接使用
2. **openpyxl**: Excel 输出 (预留)
