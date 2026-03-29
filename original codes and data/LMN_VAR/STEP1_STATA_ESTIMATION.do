**************************************************************
* STEP1_STATA_ESTIMATION.do
*
* This STATA code calls the VAR coefficient estimation
* files used for the disaster event restrictions VAR in
* "Using Disasters to Estimate the Impact of Uncertainty"
* by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
* The file is meant to be run before the 
* MATLAB files STEP2_MATLAB_ESTIMATION.m and STEP3_GRAPHS.m
* in order to produce Figures 3-5 in the paper.
*
* Initially run on STATA/MP 17.0
*
* Questions to Stephen Terry
* stephenjamesterry@gmail.com
**************************************************************

clear

*produce baseline VAR coeff matrices
do "BASELINE/var_OLSest.do"

*produce VAR coeff matrices for looser restrictions version
do "LOOSER/var_OLSest.do"

*produce VAR coeff matrices for tighter restrictions version
do "TIGHTER/var_OLSest.do"

*produce VAR coeff matrices for no country FE version
do "NO_COUNTRY_FE/var_OLSest.do"

*produce VAR coeff matrices for no time FE version
do "NO_TIME_FE/var_OLSest.do"

*produce VAR coeff matrices for 10-lag version
do "10LAGS/var_OLSest.do"

*produce VAR coeff matrices for 14-lag version
do "14LAGS/var_OLSest.do"

*produce VAR coeff matrices for revolutions only version
do "REV_ONLY/var_OLSest.do"

*produce VAR coeff matrices for revolutions-coups only version
do "REV_COUP_ONLY/var_OLSest.do"
