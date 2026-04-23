# 不会编程也能用：3步跑通本项目

这个项目用于复现论文 *Using Disasters to Estimate the Impact of Uncertainty* 的 Python 版本结果。

你只需要改一个文件：`config/experiment_config.json`。

## 1. 第一次使用（安装并测试）

在项目目录打开终端后，依次运行：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python run_all.py
```

成功后结果会写入：

- `output/runs/<任务名>/tables/`（IV 表格）
- `output/runs/<任务名>/figures/`（IV_VAR / LMN_VAR / MODEL 图形）
- `output/runs/run_summary.json`（是否成功）

## 2. 只改这一个配置文件

文件：`config/experiment_config.json`

你通常只要改以下几处：

1. `jobs`：任务列表（可以配置多组数据并行跑）
2. 每个模块的 `data_path`：输入数据路径
3. `default_modules`：要跑哪些模块（`iv` / `iv_var` / `lmn_var` / `model`）
4. `parallel.enabled` 和 `parallel.max_workers`：是否并行

改完后运行：

```bash
python run_all.py
```

只跑某些模块：

```bash
python run_all.py iv iv_var
```

只跑某个任务：

```bash
python run_all.py --jobs baseline
```

### 2.1 并行计算三种场景（可直接复制）

并行总开关（先确认）：

```json
"parallel": {
  "enabled": true,
  "max_workers": 3
}
```

说明：
- 并行是按 `jobs` 并行。
- 同时有多个 `enabled: true` 的 job 时，`python run_all.py` 会并行执行它们。

场景 A：同一数据，不同参数

```json
"jobs": [
  {
    "name": "same_data_param_a",
    "enabled": true,
    "iv_var": {
      "data_path": "data/IV_VAR/VARdata.csv",
      "bootstrap_n": 150
    }
  },
  {
    "name": "same_data_param_b",
    "enabled": true,
    "iv_var": {
      "data_path": "data/IV_VAR/VARdata.csv",
      "bootstrap_n": 300
    }
  }
]
```

场景 B：同一参数，不同数据

```json
"jobs": [
  {
    "name": "same_param_data_a",
    "enabled": true,
    "iv_var": {
      "data_path": "data/my_data/vardata_a.csv",
      "bootstrap_n": 150
    }
  },
  {
    "name": "same_param_data_b",
    "enabled": true,
    "iv_var": {
      "data_path": "data/my_data/vardata_b.csv",
      "bootstrap_n": 150
    }
  }
]
```

场景 C：不同数据，不同参数

```json
"jobs": [
  {
    "name": "diff_both_a",
    "enabled": true,
    "iv_var": {
      "data_path": "data/my_data/vardata_a.csv",
      "bootstrap_n": 100
    },
    "lmn_var": {
      "data_path": "data/my_data/lmn_a.dta",
      "n_draws": 300000
    }
  },
  {
    "name": "diff_both_b",
    "enabled": true,
    "iv_var": {
      "data_path": "data/my_data/vardata_b.csv",
      "bootstrap_n": 250
    },
    "lmn_var": {
      "data_path": "data/my_data/lmn_b.dta",
      "n_draws": 800000
    }
  }
]
```

并行计算快速示例：

```bash
# 按 jobs 并行运行
python run_all.py

# 只并行跑指定两个任务
python run_all.py --jobs same_data_param_a same_data_param_b

# 用同一配置强制串行（对照测试）
python run_all.py --sequential
```

## 3. 用项目自带数据直接出结果

默认配置已经指向项目数据：

- `data/IV/panel_iv_data.dta`
- `data/IV_VAR/VARdata.csv`
- `data/LMN_VAR/Dates_and_Data.dta`

你不改数据路径也可以直接运行并得到论文对应模块结果。

## 4. 用你自己的数据（格式要求 + 模板）

### 4.1 IV_VAR（CSV）

- 必须是 `.csv`
- 列名和顺序按模板：`data/templates/IV_VAR_template.csv`
- 关键列（共9列）：
  - `ydgdp, cs_index_ret, cs_index_vol, savgnatshock, savgpolshock, savgrevshock, savgtershock, countrycode, timecode`

### 4.2 LMN_VAR（Stata DTA）

- 必须是 `.dta`
- 模板：`data/templates/LMN_VAR_template.dta`
- 必需列：
  - `country, year, quarter, ydgdp, vol, ret`
  - 事件列至少应包含：`polshock, tershock, natshock, revshock`

### 4.3 IV（Stata DTA）

- 必须是 `.dta`
- IV 模块列很多，建议直接按完整模板改：`data/templates/IV_full_template.dta`
- 可查看 CSV 版列结构（便于人工查看列名）：`data/templates/IV_full_template.csv`

### 4.4 把你的数据接入运行

在 `config/experiment_config.json` 中，把对应模块 `data_path` 改成你的文件路径，例如：

```json
"iv_var": {
  "data_path": "data/my_data/my_iv_var.csv"
}
```

## 5. 常见报错处理

1. `python: command not found`
- 用 `python3` 代替 `python`。

2. `ModuleNotFoundError: ...`
- 先激活环境：`source .venv/bin/activate`
- 再安装依赖：`pip install -r requirements.txt`

3. `FileNotFoundError: Config not found` 或数据文件找不到
- 检查你是否在项目根目录运行。
- 检查 `experiment_config.json` 里的 `data_path` 拼写是否正确。

4. `No enabled jobs matched selection`
- 你指定了 `--jobs`，但配置里的任务 `enabled` 是 `false`，或任务名不一致。

5. `KeyError` / 列名缺失 / 读取 dta 或 csv 失败
- 你的输入数据列名不符合模板。
- 先用 `data/templates/` 下模板对齐列名和顺序，再运行。

## 6. 给完全没用过 GitHub 的用户

你可以完全不使用 GitHub 命令：

1. 直接下载本项目压缩包并解压
2. 按上面的“第1步”安装
3. 只修改 `config/experiment_config.json`
4. 运行 `python run_all.py`

---

如果你需要参数与论文逐条对应说明，请看：`config/PARAMETER_REFERENCE.md`。
