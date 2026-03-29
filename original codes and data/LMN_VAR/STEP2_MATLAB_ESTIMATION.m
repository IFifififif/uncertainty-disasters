%%%%%%%
% STEP2_MATLAB_ESTIMATION.m
%
% This MATLAB code computes admissible sets
% used for the disaster event restrictions VAR in
% "Using Disasters to Estimate the Impact of Uncertainty"
% by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
% The file is meant to be run after the 
% STATA file1 STEP1_STATA_ESTIMATION.do and before the
% MATLAB file STEP3_GRAPHS.m in order to produce Figures 3-5 in the paper.
%
% Initially run on MATLAB R2021a
%
% Questions to Stephen Terry
% stephenjamesterry@gmail.com
%%%%%%%

close all; clear all; clc;

%produce baseline VAR results
run BASELINE/var_LMN.m

%produce looser restriction VAR results
run LOOSER/var_LMN.m

%produce tighter restriction VAR results
run TIGHTER/var_LMN.m

%produce no country FE VAR results
run NO_COUNTRY_FE/var_LMN.m

%produce no time FE VAR results
run NO_TIME_FE/var_LMN.m

%produce 10-lag VAR results
run 10LAGS/var_LMN.m

%produce 14-lag VAR results
run 14LAGS/var_LMN.m

%produce revolution-only VAR results
run REV_ONLY/var_LMN.m

%produce revolution-coup-only VAR results
run REV_COUP_ONLY/var_LMN.m