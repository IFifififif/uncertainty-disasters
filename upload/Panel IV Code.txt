clear all
use panel_iv_data


****************************************************************************************************************************
****************************************************************************************************************************
*****TABLE 1: D-STATS for sample with non-missing interpolation yearly GDP growth (largest possible sample)*****************
eststo summstats: estpost summarize ydgdp cs_index_ret cs_index_vol avgret lavgvol avgcs_ret lavgcs_vol ynatshock ysavgnatshock ypolshock ysavgpolshock yrevshock ysavgrevshock ytershock ysavgtershock GDP if ydgdp~=.,de
esttab summstats using "dstats.csv", replace cell("count mean p50 sd min max")

***********************************************************************************************************************************
*****TABLE 2: BASELINE
global iv "l1savgnatshock l1savgpolshock l1savgrevshock l1savgtershock"
global d_iv "l1savgd_natshock l1savgd_polshock l1savgd_revshock l1savgd_tershock"
global t_iv "l1savgt_natshock l1savgt_polshock l1savgt_revshock l1savgt_tershock"

*Col 1 annual OLS - micro+macro
areg ydgdp cs_index_ret cs_index_vol i.yq_int ,ab(country) cluster(country)
gen sample1 = e(sample)
outreg2 using "Baseline.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES) drop(i.yq_int) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label replace

*Col 2: - yearly, IV - micro+macro
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=$iv), cluster(country) partial(yy* cc*) first savefirst
outreg2 using "Baseline.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)', F-Stat, `e(F)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 
est restore _ivreg2_cs_index_ret
outreg2 using "Baseline.xls", ctitle("ret_1st") addtext(" - Period FE", YES," - Country FE", YES) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 
est restore _ivreg2_cs_index_vol
outreg2 using "Baseline.xls", ctitle("vol_1st") addtext(" - Period FE", YES," - Country FE", YES) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 3: - yearly, IV, Stock Index
ivreg2 ydgdp cc* yy* (l1avgret l1lavgvol = $iv) , cluster(country) partial(yy* cc*) first savefirst
outreg2 using "Baseline.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)', F-Stat, `e(F)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 
est restore _ivreg2_l1avgret
outreg2 using "Baseline.xls", ctitle("ret_1st") addtext(" - Period FE", YES," - Country FE", YES) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 
est restore _ivreg2_l1lavgvol
outreg2 using "Baseline.xls", ctitle("vol_1st") addtext(" - Period FE", YES," - Country FE", YES) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 4: - yearly, IV, Stock Index - common sample
ivreg2 ydgdp cc* yy* (l1avgret l1lavgvol = $iv) if sample1==1, cluster(country) partial(yy* cc*) first savefirst
outreg2 using "Baseline.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)', F-Stat, `e(F)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 
est restore _ivreg2_l1avgret
outreg2 using "Baseline.xls", ctitle("ret_1st") addtext(" - Period FE", YES," - Country FE", YES) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 
est restore _ivreg2_l1lavgvol
outreg2 using "Baseline.xls", ctitle("vol_1st") addtext(" - Period FE", YES," - Country FE", YES) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 5: - yearly, IV, Cross Sect
ivreg2 ydgdp cc* yy* (l1avgcs_ret l1lavgcs_vol = $iv), cluster(country)  partial(yy* cc*) first savefirst
outreg2 using "Baseline.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)', F-Stat, `e(F)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 
est restore _ivreg2_l1avgcs_ret
outreg2 using "Baseline.xls", ctitle("ret_1st") addtext(" - Period FE", YES," - Country FE", YES) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 
est restore _ivreg2_l1lavgcs_vol
outreg2 using "Baseline.xls", ctitle("vol_1st") addtext(" - Period FE", YES," - Country FE", YES) alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 


************************************************************************************************************************************************************************
*****TABLE 3: ROBUSTNESS & HIGHER MOMENTS 
*Col 1: Baseline - IV, Yearly, PCF
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=$iv) , cluster(country)  partial(yy* cc*)
outreg2 using "Robustness.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label replace

*Col 2: Population weighted
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol= $iv) [aw=lpop], cluster(country) partial(yy* cc*) 
outreg2 using "Robustness.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 3: IV, Yearly, PCF - Add skewness
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol cs_index_skew =$iv) ,  partial(yy* cc*)
outreg2 using "Robustness.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 4: IV, Yearly, PCF - Keep vol and skewness
ivreg2 ydgdp cc* yy* (cs_index_vol cs_index_skew =$iv) ,  partial(yy* cc*)
outreg2 using "Robustness.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 5: HAR adjusted micro+macro vol and ret (at daily level)
ivreg2 ydgdp cc* yy* (cs_index_ret_har cs_index_vol_har=$iv), cluster(country)  partial(yy* cc*)
outreg2 using "Robustness.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 6: Residualized adjusted micro+macro vol and ret (at quarterly level)
ivreg2 ydgdp cc* yy* (cs_index_ret_har_q cs_index_vol_har_q=$iv), cluster(country)  partial(yy* cc*)
outreg2 using "Robustness.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 


************************************************************************************************************************************************************************
****Table 4 Trade and Distance Weighted Instruments
*Col 1: - yearly, IV, trade-weighted-shocks
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=$t_iv) if sample1==1, cluster(country) first partial(yy* cc*)
outreg2 using "Weighting.xls", ctitle("trade") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label replace

*Col 2: - yearly, IV, distance-weighted-shocks
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=$d_iv) if sample1==1, cluster(country) first partial(yy* cc*)
outreg2 using "Weighting.xls", ctitle("dist") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 3: - yearly, IV, trade-weighted-shocks
ivreg2 ydgdp cc* yy* (l1avgret l1lavgvol=$t_iv) if sample1==1, cluster(country) first partial(yy* cc*)
outreg2 using "Weighting.xls", ctitle("trade") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 4: - yearly, IV, distance-weighted-shocks
ivreg2 ydgdp cc* yy* (l1avgret l1lavgvol=$d_iv) if sample1==1, cluster(country) first partial(yy* cc*)
outreg2 using "Weighting.xls", ctitle("dist") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 5: - yearly, IV, trade-weighted-shocks
ivreg2 ydgdp cc* yy* (l1avgcs_ret l1lavgcs_vol=$t_iv) if sample1==1, cluster(country) first partial(yy* cc*)
outreg2 using "Weighting.xls", ctitle("trade") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Col 6: - yearly, IV, distance-weighted-shocks
ivreg2 ydgdp cc* yy* (l1avgcs_ret l1lavgcs_vol=$d_iv) if sample1==1, cluster(country) first partial(yy* cc*)
outreg2 using "Weighting.xls", ctitle("dist") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 


******************************************************************************************************************************************************************
*****TABLE 5: Unweighted by Media Citations and Jump Definitions***********************************************************************************************************************
*Baseline
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=$iv), cluster(country) partial(yy* cc*) savefirst
outreg2 using "MediaWeightings.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label replace

*Unscaled
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=l1s0avgnatshock l1s0avgpolshock l1s0avgrevshock l1s0avgtershock), cluster(country)  partial(yy* cc*)
outreg2 using "MediaWeightings.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Over median jump (1.4)
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=l1sMedavgnatshock l1sMedavgpolshock l1sMedavgrevshock l1sMedavgtershock), cluster(country)  partial(yy* cc*)
outreg2 using "MediaWeightings.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Narrower Window - 5 pre 5 post
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=l1s2avgnatshock l1s2avgpolshock l1s2avgrevshock l1s2avgtershock), cluster(country)  partial(yy* cc*)
outreg2 using "MediaWeightings.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Do jump based on regression on time trend + post dummy
*Scaled jump reg
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=l1sprdavgnatshock l1sprdavgpolshock l1sprdavgrevshock l1sprdavgtershock), cluster(country)  partial(yy* cc*)
outreg2 using "MediaWeightings.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Nonwestern Jump (10 day)
ivreg2 ydgdp cc* yy* (cs_index_ret cs_index_vol=l1nwsavgnatshock l1nwsavgpolshock l1nwsavgrevshock l1nwsavgtershock), cluster(country)  partial(yy* cc*)
outreg2 using "MediaWeightings.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 


******************************************************************************************************************************************************************
*****TABLE 6: Other Uncertainty Proxies
global iv "l1savgnatshock l1savgpolshock l1savgrevshock l1savgtershock"

*Column 1: Use WUI
ivreg2 ydgdp cc* yy* (cs_index_ret l1lavgWUI  = $iv), partial(yy* cc*) 
outreg2 using "AlternativeUncertainty.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label replace

*Column 2: Use EPU 
ivreg2 ydgdp cc* yy* (cs_index_ret l1lavgEPU  = $iv), partial(yy* cc*)
outreg2 using "AlternativeUncertainty.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label

*Column 3: Use EPU + WUI
egen any_epu = mean(EPU), by(country)
gen epu_wui = l1lavgEPU
replace epu_wui = l1lavgWUI if l1lavgEPU==. & any_epu==.
ivreg2 ydgdp cc* yy* (cs_index_ret epu_wui  = $iv), partial(yy* cc*)
outreg2 using "AlternativeUncertainty.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label

*Column 4: Use Consensus Forecasts
ivreg2 ydgdp cc* yy* (cs_index_ret l1lgdp_for_sd = $iv) ,  partial(yy* cc*)
outreg2 using "AlternativeUncertainty.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

*Column 5: Use Exchange Rates
ivreg2 ydgdp cc* yy* (cs_index_ret l1lavgexchgvol =$iv) , cluster(country) partial(yy* cc*)
outreg2 using "AlternativeUncertainty.xls", ctitle("ydgdp") addtext(" - Period FE", YES," - Country FE", YES, Hansen J, `e(jp)') alpha(0.01, 0.05, 0.1) symbol(***,**,*) bdec(3) se nocons label 

