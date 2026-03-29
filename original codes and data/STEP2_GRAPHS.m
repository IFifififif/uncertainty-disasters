%%%%%%%
% STEP2_GRAPHS.m
%
% This MATLAB code produces the IV-VAR figures contained in
% "Using Disasters to Estimate the Impact of Uncertainty,"
% by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
% The file must be run within the IV_VAR folder of the replication 
% packet,  after the file "STEP1_ESTIMATION.m," to reproduce Figures 6-7 
% in the paper.
%
% Initially run on MATLAB R2021a
%
% Questions to Stephen Terry
% stephenjamesterry@gmail.com
%%%%%%%

close all; clear all; clc;

lwidnum = 2.5;
fsizenum=14;

%read in baseline point est
load BASELINE/IV_VAR.mat
IRF_BASE = IRF_S_TO_Y;

%read in  SEs
load BOOT/IV_VAR_BOOT.mat
IRF_SE = IRF_S_TO_Y_SE;

%read in early sample est
load EARLY/IV_VAR.mat
IRF_EARLY = IRF_S_TO_Y;

%read in late sample est
load LATE/IV_VAR.mat
IRF_LATE = IRF_S_TO_Y;

%read in fewer lags est
load FEWER_LAGS/IV_VAR.mat
IRF_FEWER_LAGS = IRF_S_TO_Y;

%read in more lags est
load MORE_LAGS/IV_VAR.mat
IRF_MORE_LAGS = IRF_S_TO_Y;

%read in no country FE est
load NO_COUNTRY_FE/IV_VAR.mat
IRF_NO_COUNTRY_FE = IRF_S_TO_Y;

%read in no time FE est
load NO_TIME_FE/IV_VAR.mat
NO_TIME_FE = IRF_S_TO_Y;


%Figure 6
IRFsamp = 1:15;
figure; 
plot(IRFsamp,IRF_BASE,'bo-',...
    IRFsamp,IRF_BASE-1.645*IRF_SE,'b--',...
    IRFsamp,IRF_BASE+1.645*IRF_SE,'b--',...
    IRFsamp,0*IRFsamp,'k--','LineWidth',lwidnum)
axis([-Inf Inf -10 5])
set(gca,'FontSize',fsizenum)
xlabel('Quarters, Shock in Period 1')
ylabel('GDP Growth, Percent Year-on-Year')

saveas(gcf,'FIGURE6.pdf')

%Figure 7
figure; 
p = plot(IRFsamp,IRF_BASE-1.645*IRF_SE,'b--',...
    IRFsamp,IRF_BASE+1.645*IRF_SE,'b--',...
    IRFsamp,IRF_FEWER_LAGS,'c+-',...
    IRFsamp,IRF_MORE_LAGS,'gs-',...
    IRFsamp,NO_TIME_FE,'m>-',...
    IRFsamp,IRF_NO_COUNTRY_FE,'rx-',...
    IRFsamp,IRF_EARLY,'kh-',...
    IRFsamp,IRF_LATE,'yp-',...
    IRFsamp,IRF_BASE,'bo-',...
    IRFsamp,0*IRFsamp,'k--','LineWidth',lwidnum);
axis([-Inf Inf -10 5])
set(gca,'FontSize',fsizenum)

xlabel('Quarters, Shock in Period 1')
ylabel('GDP Growth, Percent Year-on-Year')

p(7).Color = [1 0.5 0]; %early is orange
p(6).Color = [171 104 87]./255; %no cfe is brown

saveas(gcf,'FIGURE7.pdf')
