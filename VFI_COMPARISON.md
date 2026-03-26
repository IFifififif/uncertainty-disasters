# VFI 模块 Fortran vs Python 详细对比

## 1. 数组维度

| 数组 | Fortran 维度 | Python 维度 | 一致性 |
|------|-------------|-------------|--------|
| V | (numendog, numexog, numfcst) | (numendog, numexog, kbarnum) | ✅ |
| EVmat | (numexog, numendog, numfcst) | (numexog, numendog, kbarnum) | ✅ |
| polmat | (numendog, numexog, numfcst) | (numendog, numexog, kbarnum) | ✅ |

## 2. Howard 加速步骤

### Fortran (行 814-867)
```fortran
do accelct=1,accelmaxit
    do ct=1,numstates
        ! 提取状态索引
        endogct=loopind(ct,1); exogct=loopind(ct,2); fcstct=loopind(ct,3)
        
        ! 获取当前策略
        polstar = polmatold(endogct,exogct,fcstct)
        
        ! 计算当期回报
        V(endogct,exogct,fcstct) = pfcstmat(act,sct,smin1ct,fcstct) * (...)
        
        ! 插值计算期望延续值
        Vnextval = weight * Vold(polstar,exogprimect,ind+1) + &
                   (1.0 - weight)*Vold(polstar,exogprimect,ind)
        
        V(endogct,exogct,fcstct) = V(endogct,exogct,fcstct) + &
            beta * pr_mat(exogct,exogprimect) * Vnextval
    end do
    Vold = V
end do
```

### Python
```python
for accel in range(n_accel):
    for ct in prange(numstates):
        # 获取策略
        polstar = polmat[endogct, exogct, fcstct]
        
        # 计算当期回报
        period_return = pfcstmat[act, sct, smin1ct, fcstct] * (...)
        
        # 插值计算期望延续值
        Vnextval = weight * V_old[polstar, exogprimect, ind + 1] + \
                   (1.0 - weight) * V_old[polstar, exogprimect, ind]
        
        V_new[endogct, exogct, fcstct] = period_return + beta * EV
```

**差异分析：** ✅ 逻辑一致

## 3. EVmat 创建

### Fortran (行 870-902)
```fortran
do ct=1,numstates
    polct=loopind(ct,1); exogct=loopind(ct,2); fcstct=loopind(ct,3)
    
    EVmat(exogct,polct,fcstct) = 0.0
    do exogprimect=1,numexog
        ! 关键：这里使用 polct (策略索引)，不是 polstar
        Vnextval = weight * Vold(polct,exogprimect,ind+1) + &
                   (1.0 - weight)*Vold(polct,exogprimect,ind)
        
        EVmat(exogct,polct,fcstct) = EVmat(exogct,polct,fcstct) + &
            pr_mat(exogct,exogprimect) * Vnextval
    end do
end do
```

### Python
```python
for ct in prange(numendog * numexog * kbarnum):
    polct = ct // (numexog * kbarnum)
    exogct = (ct % (numexog * kbarnum)) // kbarnum
    fcstct = ct % kbarnum
    
    for exogprimect in range(numexog_next):
        Vnextval = weight * V[polct, exogprimect, ind + 1] + \
                   (1.0 - weight) * V[polct, exogprimect, ind]
        ev += pr_mat[exogct, exogprimect] * Vnextval
    
    EVmat[exogct, polct, fcstct] = ev
```

**差异分析：** ✅ 逻辑一致

## 4. 优化步骤

### Fortran (行 904-943)
```fortran
do ct=1,numstates
    endogct=loopind(ct,1); exogct=loopind(ct,2); fcstct=loopind(ct,3)
    
    ! 遍历所有策略
    do polct=1,numendog
        ! 计算当期回报
        RHSvec(polct) = pfcstmat(act,sct,smin1ct,fcstct) * (...)
        
        ! 加上期望延续值
        RHSvec(polct) = RHSvec(polct) + beta * EVmat(exogct,polct,fcstct)
    end do
    
    ! 选择最优策略
    polstar = maxloc(RHSvec,1)
    polmat(endogct,exogct,fcstct) = polstar
    V(endogct,exogct,fcstct) = RHSvec(polstar)
end do
```

### Python
```python
for ct in prange(numstates):
    # 遍历所有策略
    for polct in range(numendog):
        period_return = pfcstmat[act, sct, smin1ct, fcstct] * (...)
        val = period_return + beta * EVmat[exogct, polct, fcstct]
        
        if val > best_val:
            best_val = val
            best_pol = polct
    
    V_new[endogct, exogct, fcstct] = best_val
    polmat_new[endogct, exogct, fcstct] = best_pol
```

**差异分析：** ✅ 逻辑一致

## 5. 索引差异处理

| 项目 | Fortran | Python | 处理方式 |
|------|---------|--------|----------|
| 数组索引 | 1-based | 0-based | 所有索引减1 |
| kbarfcstinds 初始化 | 1 | 0 | 已修复 |
| maxloc 返回值 | 1 到 n | 0 到 n-1 | argmax 返回 0-based |

## 6. 关键修复记录

### 修复 1: kbarfcstinds 初始化
- **问题**: Fortran 初始化为 1，Python 初始化为 ones 导致索引越界
- **修复**: 改为 zeros 初始化

### 修复 2: 边界检查
- **问题**: ind+1 可能超出 kbarnum 范围
- **修复**: 添加边界检查 `if ind >= kbarnum - 1`

## 7. 测试结果

| 测试项 | 结果 |
|--------|------|
| VFI 收敛 | ✅ 12 次迭代收敛 |
| VF 误差 | ✅ ~1e-12 |
| 策略误差 | ✅ 0 (收敛) |
| 边界警告 | ⚠️ 资本/劳动网格边界 |

