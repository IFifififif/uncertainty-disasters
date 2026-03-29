%%%%%%%
% STEP3_GRAPHS.m
%
% This MATLAB code compiles results into figures for 
% the disaster event restrictions VAR in
% "Using Disasters to Estimate the Impact of Uncertainty"
% by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
% The file is meant to be run after the 
% STATA file1 STEP1_STATA_ESTIMATION.do and the
% MATLAB file STEP2_MATLAB_ESTIMATION.m in order to produce Figures 3-5 in the paper.
%
% Initially run on MATLAB R2021a
%
% Questions to Stephen Terry
% stephenjamesterry@gmail.com
%%%%%%%

close all; clear all; clc

lwidnum = 0.5; %line width
fsizenum = 14; %font size
Tirf = 50; %number of periods to plot

%load baseline VAR results
load BASELINE/LMN_WORK.mat
IRF_admit_lb_BASE = IRF_admit_lb;
IRF_admit_ub_BASE = IRF_admit_ub;
IRF_maxg_BASE = IRF_maxg;
IRF_med_BASE = IRF_med;
Impact_Hist_BASE = IMPACT_HIST;

%load revolutions-coups only VAR results
load REV_COUP_ONLY/LMN_WORK.mat
IRF_admit_lb_REMOVE_NAT_TER = IRF_admit_lb;
IRF_admit_ub_REMOVE_NAT_TER = IRF_admit_ub;
Impact_Hist_REMOVE_NAT_TER = IMPACT_HIST;

%load revolutions-only VAR results
load REV_ONLY/LMN_WORK.mat
IRF_admit_lb_REMOVE_NAT_TER_COUP = IRF_admit_lb;
IRF_admit_ub_REMOVE_NAT_TER_COUP = IRF_admit_ub;
Impact_Hist_REMOVE_NAT_TER_COUP = IMPACT_HIST;

%load tighter restrictions VAR results
load TIGHTER/LMN_WORK.mat
IRF_admit_lb_TIGHTER = IRF_admit_lb;
IRF_admit_ub_TIGHTER = IRF_admit_ub;

%load looser restrictions VAR results
load LOOSER/LMN_WORK.mat
IRF_admit_lb_LOOSER = IRF_admit_lb;
IRF_admit_ub_LOOSER = IRF_admit_ub;

%load 10-lag VAR results
load 10LAGS/LMN_WORK.mat
IRF_admit_lb_10LAGS = IRF_admit_lb;
IRF_admit_ub_10LAGS = IRF_admit_ub;

%load 14-lag VAR results
load 14LAGS/LMN_WORK.mat
IRF_admit_lb_14LAGS = IRF_admit_lb;
IRF_admit_ub_14LAGS = IRF_admit_ub;

%load no country FE VAR results
load NO_COUNTRY_FE/LMN_WORK.mat
IRF_admit_lb_NOCFE = IRF_admit_lb;
IRF_admit_ub_NOCFE = IRF_admit_ub;

%load no time FE VAR results
load NO_TIME_FE/LMN_WORK.mat
IRF_admit_lb_NOTFE = IRF_admit_lb;
IRF_admit_ub_NOTFE = IRF_admit_ub;

IRFaxis = [ 1 15 -4 1];

%FIGURE 3 - BASELINE
figure;
plot(1:Tirf,squeeze(IRF_admit_lb_BASE(1,3,:)),'b',...
1:Tirf,squeeze(IRF_admit_ub_BASE(1,3,:)),'b',...
1:Tirf,squeeze(IRF_maxg_BASE(1,3,:)),'r-o',...
1:Tirf,squeeze(IRF_med_BASE(1,3,:)),'g-x',...
1:Tirf,0*(1:Tirf),'k--',...
    'LineWidth',lwidnum);
axis(IRFaxis)
set(gca,'FontSize',fsizenum);
xlabel('Quarters, Shock in Period 1');
ylabel('GDP Growth, Percent Year-on-Year')
saveas(gcf,'FIGURE3.pdf')


%FIGURE 4 - IMPACT HISTOGRAMS
figure;
h4=histogram(Impact_Hist_REMOVE_NAT_TER_COUP,25); hold on;
h4.Normalization = 'probability'; h4.FaceColor=[1 0 0];  h4.FaceAlpha=0.8;

h2=histogram(Impact_Hist_REMOVE_NAT_TER,25); 
h2.Normalization = 'probability'; h2.FaceColor=[0 1 0];  h2.FaceAlpha=0.8;

h1=histogram(Impact_Hist_BASE,25);
h1.Normalization = 'probability'; h1.FaceColor=[0 0 1]; h1.FaceAlpha=0.8;
xlabel('GDP Growth Response to Uncertainty');
ylabel('Probability')
set(gca,'FontSize',fsizenum);
saveas(gcf,'FIGURE4.pdf')

%FIGURE 5 - ROBUSTNESS
figure;
plot(1:Tirf,squeeze(IRF_admit_lb_BASE(1,3,:)),'b',...
1:Tirf,squeeze(IRF_admit_ub_BASE(1,3,:)),'b',...
1:Tirf,squeeze(IRF_maxg_BASE(1,3,:)),'r-o',...
1:Tirf,squeeze(IRF_admit_lb_14LAGS(1,3,:)),'g-x',...
1:Tirf,squeeze(IRF_admit_ub_14LAGS(1,3,:)),'g-x',...
1:Tirf,squeeze(IRF_admit_lb_TIGHTER(1,3,:)),'m-d',...
1:Tirf,squeeze(IRF_admit_ub_TIGHTER(1,3,:)),'m-d',...
1:Tirf,squeeze(IRF_admit_lb_LOOSER(1,3,:)),'c-*',...
1:Tirf,squeeze(IRF_admit_ub_LOOSER(1,3,:)),'c-*',...
1:Tirf,squeeze(IRF_admit_lb_NOTFE(1,3,:)),'y-s',...
1:Tirf,squeeze(IRF_admit_ub_NOTFE(1,3,:)),'y-s',...
1:Tirf,0*(1:Tirf),'k--',...
    'LineWidth',lwidnum); hold on;
axis(IRFaxis)
plot(1:Tirf,squeeze(IRF_admit_lb_NOCFE(1,3,:)),'Color',[0.4940 0.1840 0.5560],'Marker','^','LineWidth',lwidnum);
plot(1:Tirf,squeeze(IRF_admit_ub_NOCFE(1,3,:)),'Color',[0.4940 0.1840 0.5560],'Marker','^','LineWidth',lwidnum);
plot(1:Tirf,squeeze(IRF_admit_lb_10LAGS(1,3,:)),'Color',[0.91 0.41 0.14],'Marker','o','LineWidth',lwidnum);
plot(1:Tirf,squeeze(IRF_admit_ub_10LAGS(1,3,:)),'Color',[0.91 0.41 0.14],'Marker','o','LineWidth',lwidnum);

set(gca,'FontSize',fsizenum);
xlabel('Quarters, Shock in Period 1');
ylabel('GDP Growth, Percent Year-on-Year')
saveas(gcf,'FIGURE5.pdf')


%FIGURE B1 (NOT IN PUBLICATION)
shocklabelsvec = {'Growth Shock','First Shock','Second Shock'};
labelsvec = {'Growth','First','Second'};
IRFaxis = [ 1 15 -Inf Inf];

figure;
plotct=0;
for shockct=1:3    
    for varct=1:3
        plotct=plotct+1;
        subplot(3,3,plotct)
        plot(1:Tirf,0*(1:Tirf),'k--',...
        1:Tirf,squeeze(IRF_admit_lb_BASE(varct,shockct,:)),'b',...
        1:Tirf,squeeze(IRF_admit_ub_BASE(varct,shockct,:)),'b',...
        'LineWidth',lwidnum);
    axis(IRFaxis)

    if (shockct==1)
       title(labelsvec{varct})
    end
    if (shockct==3)
        xlabel('Quarters');
    end
    if (varct==1)
        ylabel(shocklabelsvec{shockct});
    end

    set(gca,'FontSize',fsizenum);

    end
end
saveas(gcf,'FIGUREB1.pdf')


