# Quick Check Discrepancy Table

| module   | metric                         | reference                            |   python_value | delta                  | status   |
|:---------|:-------------------------------|:-------------------------------------|---------------:|:-----------------------|:---------|
| IV       | Table2 Col2 coef(cs_index_ret) | Stata ivreg2 output                  |      1.19684   | N/A (no baseline file) | run-ok   |
| IV       | Table2 Col2 coef(cs_index_vol) | Stata ivreg2 output                  |     -4.23631   | N/A (no baseline file) | run-ok   |
| IV_VAR   | Baseline impact t=1            | MATLAB VAR.m (negative GDP response) |     -0.0963931 | sign-ok                | run-ok   |
| IV_VAR   | Bootstrap SE finite check      | finite SE vector                     |      9.72565   | ok                     | run-ok   |
| LMN_VAR  | Execution status               | must run end-to-end                  |      1         | ok                     | run-ok   |
| LMN_VAR  | Admissible share (5k draws)    | >0 often requires large N draws      |      0.0038    | ok                     | run-ok   |
| MODEL    | Moment RMSE (20 moments)       | 0 (perfect match target)             |      4.92058   | 4.920584150735211      | run-ok   |
| MODEL    | Moment MAE (20 moments)        | 0 (perfect match target)             |      3.01952   | 3.019518343710499      | run-ok   |
