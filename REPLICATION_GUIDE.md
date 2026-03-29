# 复现操作指南

## Using Disasters to Estimate the Impact of Uncertainty
### Baker, Bloom, and Terry (2024), Review of Economic Studies

---

## 目录

1. [环境要求](#1-环境要求)
2. [项目克隆](#2-项目克隆)
3. [依赖安装](#3-依赖安装)
4. [数据文件](#4-数据文件)
5. [运行IV模块 (Tables 1-6)](#5-运行iv模块-tables-1-6)
6. [运行IV_VAR模块 (Figures 6-7)](#6-运行iv_var模块-figures-6-7)
7. [运行LMN_VAR模块 (Figures 3-5)](#7-运行lmn_var模块-figures-3-5)
8. [运行MODEL模块 (Figure 8)](#8-运行model模块-figure-8)
9. [一键运行全部](#9-一键运行全部)
10. [预期输出](#10-预期输出)
11. [常见问题](#11-常见问题)

---

## 1. 环境要求

### 系统要求
- **操作系统**: Linux / macOS / Windows (WSL2)
- **Python**: 3.9 - 3.11
- **内存**: ≥8GB RAM (推荐 16GB)
- **磁盘空间**: ≥2GB

### 必需软件
```bash
# 检查 Python 版本
python --version  # 应显示 3.9.x 或更高

# 检查 pip
pip --version
```

---

## 2. 项目克隆

### 方法一：HTTPS 克隆（推荐）

```bash
# 进入工作目录
cd ~/projects

# 克隆仓库
git clone https://github.com/IFifififif/uncertainty-disasters.git

# 进入项目目录
cd uncertainty-disasters
```

### 方法二：SSH 克隆

```bash
git clone git@github.com:IFifififif/uncertainty-disasters.git
cd uncertainty-disasters
```

### 方法三：下载 ZIP

```bash
# 下载并解压
wget https://github.com/IFifififif/uncertainty-disasters/archive/refs/heads/main.zip
unzip main.zip
cd uncertainty-disasters-main
```

---

## 3. 依赖安装

### 步骤 3.1：创建虚拟环境（推荐）

```bash
# 创建虚拟环境
python -m venv .venv

# 激活虚拟环境
# Linux/macOS:
source .venv/bin/activate

# Windows:
# .venv\Scripts\activate
```

### 步骤 3.2：安装依赖

```bash
# 升级 pip
pip install --upgrade pip

# 安装项目依赖
pip install -r requirements.txt
```

### 步骤 3.3：验证安装

```bash
# 验证核心包
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"
python -c "import pandas; print(f'Pandas: {pandas.__version__}')"
python -c "import scipy; print(f'SciPy: {scipy.__version__}')"
python -c "import linearmodels; print('Linearmodels: OK')"
```

**预期输出**:
```
NumPy: 1.26.x
Pandas: 2.x.x
SciPy: 1.11.x
Linearmodels: OK
```

---

## 4. 数据文件

### 数据文件结构

```
data/
├── IV/
│   ├── panel_iv_data.dta    # IV 面板数据
│   └── dstats.csv           # 描述性统计
├── IV_VAR/
│   └── VARdata.csv          # VAR 数据
├── LMN_VAR/
│   └── Dates_and_Data.dta   # LMN VAR 数据
└── MODEL/                    # 模型生成数据
```

### 验证数据完整性

```bash
# 检查数据文件是否存在
ls -la data/IV/
ls -la data/IV_VAR/
ls -la data/LMN_VAR/
```

**注意**: 数据文件已包含在仓库中，无需额外下载。

---

## 5. 运行IV模块 (Tables 1-6)

### 步骤 5.1：运行完整IV分析

```bash
# 确保在项目根目录
cd ~/projects/uncertainty-disasters

# 运行 IV 模块
python -m src.iv.panel_iv
```

### 步骤 5.2：预期运行时间

- Table 1 (描述性统计): ~5秒
- Table 2 (基准回归): ~30秒
- Table 3 (稳健性检验): ~1分钟
- Table 4 (加权工具变量): ~1分钟
- Table 5 (媒体权重): ~1分钟
- Table 6 (替代代理): ~1分钟

**总运行时间**: ~5分钟

### 步骤 5.3：输出文件

```bash
# 查看输出
ls -la output/tables/

# 预期文件:
# table1_dstats.csv
# table2_baseline.csv
# table3_robustness.csv
# table4_weighting.csv
# table5_media_weightings.csv
# table6_alternative_proxies.csv
```

### 步骤 5.4：验证结果

```bash
# 查看Table 2基准回归结果
cat output/tables/table2_baseline.csv
```

**预期输出示例**:
```
,Col1,Col2,Col3,Col4,Col5
cs_index_ret,0.0512,0.0487,...
cs_index_vol,-0.0234,-0.0198,...
...
```

---

## 6. 运行IV_VAR模块 (Figures 6-7)

### 步骤 6.1：运行IV-VAR估计

```bash
# 运行 IV_VAR 模块
python -m src.iv_var.estimation
```

### 步骤 6.2：预期运行时间

- GMM优化: ~2-5分钟
- Bootstrap标准误 (500次): ~10-15分钟
- IRF计算: ~1分钟
- 图形生成: ~10秒

**总运行时间**: ~15-20分钟

### 步骤 6.3：输出文件

```bash
# 查看图形输出
ls -la output/figures/

# 预期文件:
# FIGURE6.pdf  # 基准 IRF
# FIGURE7.pdf  # 稳健性 IRF
```

### 步骤 6.4：验证结果

```bash
# 查看 GMM 估计结果（终端输出）
# 预期输出:
# GMM optimization converged
# Parameter estimates: [...]
# IRF computed successfully
```

---

## 7. 运行LMN_VAR模块 (Figures 3-5)

### 步骤 7.1：运行LMN VAR估计

```bash
# 运行 LMN_VAR 模块
python -m src.lmn_var.estimation
```

### 步骤 7.2：预期运行时间

- VAR FE估计: ~10秒
- 容许集计算: ~30秒
- IRF计算: ~1分钟
- 图形生成: ~10秒

**总运行时间**: ~2-3分钟

### 步骤 7.3：输出文件

```bash
# 查看图形输出
ls -la output/figures/

# 预期文件:
# FIGURE3.pdf  # GDP IRF
# FIGURE4.pdf  # 波动率 IRF
# FIGURE5.pdf  # 收益率 IRF
```

---

## 8. 运行MODEL模块 (Figure 8)

### 步骤 8.1：选择网格规模

```bash
# 完整网格 (与论文一致，需要较长时间)
# 在 src/model/params.py 中确认:
# znum = 9, anum = 21, knum = 150, lnum = 75
```

或使用简化网格（快速测试）:
```python
# 在代码中使用:
params = create_params(simplified=True)  # znum=5, anum=7, knum=30, lnum=15
```

### 步骤 8.2：运行模型求解

```bash
# 运行完整模型
python -m src.model
```

### 步骤 8.3：预期运行时间

| 步骤 | 完整网格 | 简化网格 |
|------|---------|---------|
| 网格构建 | ~5秒 | ~1秒 |
| VFI迭代 | ~10-15分钟 | ~1-2分钟 |
| 一般均衡 | ~5分钟 | ~30秒 |
| 多企业模拟 | ~2分钟 | ~10秒 |
| GMM估计 | ~10分钟 | ~1分钟 |
| IRF计算 | ~5分钟 | ~30秒 |

**总运行时间**: 
- 完整网格: ~30-40分钟
- 简化网格: ~5分钟

### 步骤 8.4：输出文件

```bash
# 查看图形输出
ls -la output/figures/

# 预期文件:
# FIGURE8.pdf  # 模型 IRF
```

### 步骤 8.5：验证结果

```bash
# 终端输出预期:
# VFI converged in 12 iterations
# GE solver completed: Y mean = 0.0098
# Multi-firm simulation: 800 firms, GDP mean = 138.10
# GMM objective value = 20652
# IRF computed: Y range [0, 1.99], I range [0, 29.16]
```

---

## 9. 一键运行全部

### 方法一：使用Python脚本

```bash
# 运行所有模块
python run_all.py
```

### 方法二：使用Shell脚本

```bash
# 创建运行脚本
cat > run_replication.sh << 'EOF'
#!/bin/bash
set -e

echo "=========================================="
echo "Starting Full Replication"
echo "=========================================="

echo ""
echo "[1/4] Running IV Module (Tables 1-6)..."
python -m src.iv.panel_iv
echo "IV Module completed."

echo ""
echo "[2/4] Running IV_VAR Module (Figures 6-7)..."
python -m src.iv_var.estimation
echo "IV_VAR Module completed."

echo ""
echo "[3/4] Running LMN_VAR Module (Figures 3-5)..."
python -m src.lmn_var.estimation
echo "LMN_VAR Module completed."

echo ""
echo "[4/4] Running MODEL Module (Figure 8)..."
python -m src.model
echo "MODEL Module completed."

echo ""
echo "=========================================="
echo "Full Replication Completed!"
echo "=========================================="
echo ""
echo "Output files:"
ls -la output/tables/
ls -la output/figures/
EOF

# 添加执行权限
chmod +x run_replication.sh

# 运行
./run_replication.sh
```

---

## 10. 预期输出

### 10.1 表格输出 (output/tables/)

| 文件 | 内容 | 论文对应 |
|------|------|---------|
| `table1_dstats.csv` | 描述性统计 | Table 1 |
| `table2_baseline.csv` | 基准IV回归 | Table 2 |
| `table3_robustness.csv` | 稳健性检验 | Table 3 |
| `table4_weighting.csv` | 贸易/距离加权 | Table 4 |
| `table5_media_weightings.csv` | 媒体权重 | Table 5 |
| `table6_alternative_proxies.csv` | 替代代理 | Table 6 |

### 10.2 图形输出 (output/figures/)

| 文件 | 内容 | 论文对应 |
|------|------|---------|
| `FIGURE3.pdf` | GDP脉冲响应 | Figure 3 |
| `FIGURE4.pdf` | 波动率脉冲响应 | Figure 4 |
| `FIGURE5.pdf` | 收益率脉冲响应 | Figure 5 |
| `FIGURE6.pdf` | IV-VAR基准IRF | Figure 6 |
| `FIGURE7.pdf` | IV-VAR稳健性IRF | Figure 7 |
| `FIGURE8.pdf` | 模型IRF | Figure 8 |

### 10.3 结果验证清单

#### Table 2 基准回归验证
```python
# 预期系数符号:
# cs_index_ret: 正 (收益率对GDP增长有正向影响)
# cs_index_vol: 负 (波动率对GDP增长有负向影响)

# 第二阶段预期结果:
# 第一矩系数: 正 (~1.5)
# 第二矩系数: 负 (~-3.8)
```

#### Figure 8 IRF验证
```python
# 预期IRF模式:
# 不确定性冲击后:
# - GDP下降 (负向响应)
# - 投资下降 (负向响应)
# - 波动率上升 (正向响应)
# 
# 响应持续约 20-30 期后回归稳态
```

---

## 11. 常见问题

### Q1: 内存不足错误

**问题**: `MemoryError` 或系统卡死

**解决方案**:
```python
# 使用简化网格
params = create_params(simplified=True)
```

### Q2: 数据文件找不到

**问题**: `FileNotFoundError: panel_iv_data.dta`

**解决方案**:
```bash
# 检查数据目录
ls data/IV/

# 如果缺失，从仓库重新获取
git checkout data/
```

### Q3: linearmodels 导入错误

**问题**: `ModuleNotFoundError: No module named 'linearmodels'`

**解决方案**:
```bash
pip install linearmodels
```

### Q4: VFI 不收敛

**问题**: `VFI did not converge after 50 iterations`

**解决方案**:
```python
# 增加迭代次数
params = ModelParameters(vfmaxit=100)

# 或放宽收敛容差
params = ModelParameters(vferrortol=1e-3)
```

### Q5: 图形显示乱码

**问题**: PDF中文字体问题

**解决方案**:
```python
# 在代码中添加:
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
```

### Q6: 并行计算警告

**问题**: `OpenMP warning` 或 `Numba警告`

**解决方案**:
```bash
# 设置环境变量
export OMP_NUM_THREADS=4
export NUMBA_NUM_THREADS=4
```

---

## 12. 进阶用法

### 12.1 修改参数进行反事实实验

```python
# 创建自定义参数
from src.model.params import ModelParameters

params = ModelParameters(
    # 修改调整成本
    capirrev=0.5,     # 增加资本不可逆成本
    hirelin=0.1,      # 增加雇佣成本
    
    # 修改不确定性持续性
    uncpers=0.9,      # 降低不确定性持续性
    
    # 修改网格大小
    knum=100,         # 减少资本网格
    lnum=50           # 减少劳动网格
)

# 运行模型
from src.model.solve import solve_model
results = solve_model(params)
```

### 12.2 仅运行特定模块

```python
# 仅运行 VFI
from src.model.vfi import run_vfi
V, policies = run_vfi(params)

# 仅运行模拟
from src.model.simulation import run_simulation
sim_results = run_simulation(params, policies)

# 仅计算 IRF
from src.model.irf import compute_irf
irf = compute_irf(params, policies)
```

### 12.3 保存和加载结果

```python
import pickle

# 保存结果
with open('results.pkl', 'wb') as f:
    pickle.dump(results, f)

# 加载结果
with open('results.pkl', 'rb') as f:
    results = pickle.load(f)
```

---

## 13. 联系与支持

### 问题报告
- GitHub Issues: https://github.com/IFifififif/uncertainty-disasters/issues

### 原始论文
- Baker, S. R., Bloom, N., & Terry, S. J. (2024). Using Disasters to Estimate the Impact of Uncertainty. *The Review of Economic Studies*, 91(2), 720–747.
- DOI: https://doi.org/10.1093/restud/rdad036

### 原始代码
- 作者提供代码: https://www.dropbox.com/s/xxx (参见论文附录)

---

**文档版本**: 1.0
**最后更新**: 2025年3月
