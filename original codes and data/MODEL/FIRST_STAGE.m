%%%%%%%
%FIRST_STAGE.m
%
% This code computes target IV regression coefficients
% estimated on simulated data as part of the structural 
% estimation of the model in 
% "Using Disasters to Estimate the Impact of Uncertainty,"
% by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry
%
% The code is called by the Fortran estimation wrapper.
%
%Questions to Stephen Terry
%stephenjamesterry@gmail.com
%%%%%%%

close all; clear all; clc;
close all; clear all; clc; 
warning('off','all')
DATA = importdata('MATLABdata.csv');
delete('FIRST_STAGE_OUT.txt');

DATAMAT = DATA.data;
%1 = placeholder var
%2 = growthsim (in %)
%3 = growthsimyr (in %)
%4 = fret macro
%5 = sret macro
%6 = fret micro
%7 = sret micro
%8 = natdissim (natural disaster shock)
%9 = polshocksim (coup shock)
%10 = revolutionsim (revolution shock)
%11 = terrorsim (terror shock)
%12 = country (categorical, not dummy)
%13 = time period

%scaled first moments by their SD
DATAMAT(:,4) = DATAMAT(:,4)/std(DATAMAT(:,4),1);
DATAMAT(:,6) = DATAMAT(:,6)/std(DATAMAT(:,6),1);

%convert second moments to logs and then scale by SD
DATAMAT(:,5) = log(DATAMAT(:,5));
DATAMAT(:,5) = DATAMAT(:,5)/std(DATAMAT(:,5),1);

DATAMAT(:,7) = log(DATAMAT(:,7));
DATAMAT(:,7) = DATAMAT(:,7)/std(DATAMAT(:,7),1);

%create dummy variables
countrydummat = dummyvar(DATAMAT(:,12));
timedummat = dummyvar(DATAMAT(:,13));

Z = [DATAMAT(:,[8 9 10 11]) countrydummat timedummat];

%macro first stage coefficients
FRETmacrobetaFIRST = (Z'*Z)\(Z'*DATAMAT(:,4)); FRETmacrobetaFIRSTdis=FRETmacrobetaFIRST(1:4);
SRETmacrobetaFIRST = (Z'*Z)\(Z'*DATAMAT(:,5)); SRETmacrobetaFIRSTdis=SRETmacrobetaFIRST(1:4);

%macro first stage predictions
FRETmacrohat = Z*FRETmacrobetaFIRST;
SRETmacrohat = Z*SRETmacrobetaFIRST;

%macro second stage coeffs
Xhat = [FRETmacrohat SRETmacrohat countrydummat timedummat];
betaSECONDmacro = (Xhat'*Xhat)\(Xhat'*DATAMAT(:,3));
betaSECONDmacro = betaSECONDmacro(1:2);

%micro first stage coefficients
FRETmicrobetaFIRST = (Z'*Z)\(Z'*DATAMAT(:,6)); FRETmicrobetaFIRSTdis=FRETmicrobetaFIRST(1:4);
SRETmicrobetaFIRST = (Z'*Z)\(Z'*DATAMAT(:,7)); SRETmicrobetaFIRSTdis=SRETmicrobetaFIRST(1:4);

%micro first stage predictions
FRETmicrohat = Z*FRETmicrobetaFIRST;
SRETmicrohat = Z*SRETmicrobetaFIRST;

%micro second stage coeffs
Xhat = [FRETmicrohat SRETmicrohat countrydummat timedummat];
betaSECONDmicro = (Xhat'*Xhat)\(Xhat'*DATAMAT(:,3));
betaSECONDmicro = betaSECONDmicro(1:2);

%micro OLS coeffs
Xhat = [DATAMAT(:,6) DATAMAT(:,7) countrydummat timedummat];
betaOLSmicro = (Xhat'*Xhat)\(Xhat'*DATAMAT(:,3));
betaOLSmicro = betaOLSmicro(1:2);

OUTVEC = [FRETmacrobetaFIRSTdis; SRETmacrobetaFIRSTdis; betaSECONDmacro; ...
   FRETmicrobetaFIRSTdis; SRETmicrobetaFIRSTdis; betaSECONDmicro];

OUTFILE_STR = ['FIRST_STAGE_OUT.txt'];
save(OUTFILE_STR,'OUTVEC','-ascii');


