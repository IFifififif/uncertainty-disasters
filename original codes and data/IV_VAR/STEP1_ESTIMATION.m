%%%%%%%
% STEP1_ESTIMATION.m
%
% This MATLAB code estimates the various IV-VAR specification in 
% "Using Disasters to Estimate the Impact of Uncertainty,"
% by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry. The file
% must be run within the IV_VAR folder of the replication packet, 
% before the file "STEP2_GRAPHS.m," to reproduce Figures 6-7 
% in the paper.
%
% Initially run on MATLAB R2021a
%
% Questions to Stephen Terry
% stephenjamesterry@gmail.com
%%%%%%%

close all; clear all; clc;

%call baseline point est
run BASELINE/VAR.m

%call boostrap to compute SEs
run BOOT/VAR.m

%call early sample est
run EARLY/VAR.m

%call late sample est
run LATE/VAR.m

%call fewer lags est
run FEWER_LAGS/VAR.m

%call more lags est
run MORE_LAGS/VAR.m

%call no country FE est
run NO_COUNTRY_FE/VAR.m

%call no time FE est
run NO_TIME_FE/VAR.m