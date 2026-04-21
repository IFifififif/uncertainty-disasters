# Paper + Original Code Gap Table (Latest)

| 模块 | 论文/原始代码目标 | 当前结果 | 结论 |
|---|---|---|---|
| IV Table 2 Col(2) | level≈1.197, vol≈-4.236（论文 Table 2） | level=1.197, vol=-4.236, Jp=0.367 | 基本一致 |
| IV_VAR Figure 6 | 不确定性冲击后 GDP 立即下滑（方向与量级） | impact_t1=-4.380 | 方向一致，量级偏深 |
| LMN_VAR Figure 3 | admissible 边界约 1%~2.5% 下滑（论文文本） | med_t1=-1.943, maxg_t1=-2.234 | 区间接近 |
| MODEL Figure 8 (smoke) | 按 Fortran 机制输出相对基线偏离 | gdp_min=-1.041, inv_min=-14.376 | 口径已对齐（需继续用 full 规模校准幅度） |
| MODEL Figure 8 (medium) | 同上 | gdp_min=-0.493, inv_min=-11.475 | 中等网格下稳定负向响应 |

## 说明
- 本轮对比基于纯 Python 流水线，不依赖 MATLAB/Fortran 文本中间文件。
- `quick_check` 默认仅输出 JSON；CSV/Markdown 为可选导出。
