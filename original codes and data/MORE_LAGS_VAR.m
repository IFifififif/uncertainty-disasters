%%%%%%%
% VAR.m
%
% This MATLAB code estimates the more lags IV-VAR in
% "Using Disasters to Estimate the Impact of Uncertainty,"
% by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
%
% Initially run on MATLAB R2021a
%
% Questions to Stephen Terry
% stephenjamesterry@gmail.com
%%%%%%%

close all; clear all; clc;
global NX ND Nparams Nmoms MOMvec extraoutput

s = RandStream('mt19937ar','Seed',3991);
RandStream.setGlobalStream(s);

Nlags = 12; 
   
fminconopt = optimoptions(@fmincon,'MaxFunEvals',50000);   
paramguess= zeros(Nparams,1);
paramguess(1) = 1; 
paramguess(2) = 0;
paramguess(3) = 0;
paramguess(4) = 0; 
paramguess(5) = 1;
paramguess(6) = 0;
paramguess(7) = 0; 
paramguess(8) = 0;
paramguess(9) = 1;
paramguess(10) = -1;
paramguess(11) = -1;
paramguess(12) = -1;
paramguess(13) = -1;
paramguess(14) = 0;
paramguess(15) = 1; 
paramguess(16) = 1;
paramguess(17) = 1;
paramguess = 0.25*paramguess;

%set basic dimensions
NX = 3; %number of variables in VAR
ND = 4; %number of instruments in IV approach
lengthIRF = 50; %length of IRFs to compute

Nparams = NX*NX + ND*2; %number of parameters to estimate
Nmoms = NX*(NX+1)/2 + ND*NX; %number of moments

%read in the data
DATACELLS = importdata('../VARdata.csv');
DATA = DATACELLS.data;

Growth = DATA(:,1); 
First = DATA(:,2); 
Second = DATA(:,3); 
NatDis = DATA(:,4);
PolShock = DATA(:,5);
Revolution = DATA(:,6);
Terror = DATA(:,7);
Country = DATA(:,8);
Time = DATA(:,9);

%create list of countries
CountryList = unique(Country);
Ncountries = length(CountryList);

X = [];
Xmin1 = [];
D = [];

CountryFlags = [];
TimeFlags = [];
for countryct=1:Ncountries;
 
    countrysamp = Country==CountryList(countryct);
    
    rawTimecountry = Time(countrysamp);
    rawCountrycountry = Country(countrysamp);
    
    rawXcountry = [Growth(countrysamp) First(countrysamp) Second(countrysamp)];
    rawDcountry =  [NatDis(countrysamp) PolShock(countrysamp) Revolution(countrysamp) Terror(countrysamp)];
    sizeXcountry = size(rawXcountry);
    Tcountry = sizeXcountry(1);
    sampleVARcountry = ((Nlags+1):Tcountry)';

    %scale first and second moments by country
    rawXcountry(sampleVARcountry,2) = rawXcountry(sampleVARcountry,2)/sqrt(var(rawXcountry(sampleVARcountry,2)));
    rawXcountry(sampleVARcountry,3) = rawXcountry(sampleVARcountry,3)/sqrt(var(rawXcountry(sampleVARcountry,3)));
    
    numobs = length(sampleVARcountry);
    samplelagsVARcountry = zeros(numobs,Nlags);
    for lagct=1:Nlags;
        samplelagsVARcountry(:,lagct) = sampleVARcountry-lagct;
    end;
    
    
    if (length(sampleVARcountry)>(Nlags+1));
    X = [X; rawXcountry(sampleVARcountry,:)];
    D = [D; rawDcountry(sampleVARcountry,:)];
  
    
    CountryFlags = [CountryFlags; rawCountrycountry(sampleVARcountry)];
    TimeFlags = [TimeFlags; rawTimecountry(sampleVARcountry)];
        
    for lagct=1:Nlags;
        if (lagct==1);
        Xmin1country = rawXcountry(samplelagsVARcountry(:,lagct),:);
        else
        Xmin1country = [Xmin1country rawXcountry(samplelagsVARcountry(:,lagct),:)];
        end;
    end;
    
    Xmin1 = [Xmin1; Xmin1country];
      end
end;

%now, compute deviations within-country
for countryct=1:Ncountries;
    countrysamp = CountryFlags==CountryList(countryct);
    numobs = sum(countrysamp);
    X(countrysamp,:) = X(countrysamp,:) - repmat(mean(X(countrysamp,:),1),numobs,1);
    Xmin1(countrysamp,:) = Xmin1(countrysamp,:) - repmat(mean(Xmin1(countrysamp,:),1),numobs,1);
    

end;

%now, compute deviations within-time
TimeList = unique(TimeFlags);
Ntimes = length(TimeList);

for ct=1:Ntimes;
    
    t = TimeList(ct);
    
    timesamp = TimeFlags==t;
    numobs = sum(timesamp);
    X(timesamp,:) = X(timesamp,:) - repmat(mean(X(timesamp,:),1),numobs,1);
    Xmin1(timesamp,:) = Xmin1(timesamp,:) - repmat(mean(Xmin1(timesamp,:),1),numobs,1);
    
    
end;

T = size(X,1);

%now, demean instruments to unit st dev and 0 mean within sample
D = D-repmat(mean(D,1),T,1);
D = D./repmat(var(D,1).^0.5,T,1);
D(isnan(D))=0;

%estimate OLS params of VAR coefficient matrices
B1hat = zeros(NX,NX*Nlags);
for varct=1:NX;
    betas = (Xmin1'*Xmin1)\(Xmin1'*X(:,varct));
    B1hat(varct,:) = betas';
end;

%compute reduced-form resids & within-sample covariance matrix
etahat = zeros(T,NX);
for t=1:T;
    etadum = X(t,:)' - B1hat*Xmin1(t,:)';
    etahat(t,:) = etadum';
end;
OMEGAhat = cov(etahat);

%compute within-sample instrument moments
EDetahat = zeros(NX,ND);
for IVct=1:ND;
for varct=1:NX;
    EDetahat(varct,IVct) = mean(etahat(:,varct).*D(:,IVct));
end;
end;

%create moment vector
MOMvec = zeros(Nmoms,1);
MOMvec(1) = OMEGAhat(1,1);
MOMvec(2) = OMEGAhat(2,2);
MOMvec(3) = OMEGAhat(3,3);
MOMvec(4) = OMEGAhat(1,2);
MOMvec(5) = OMEGAhat(1,3);
MOMvec(6) = OMEGAhat(2,3);
MOMvec(7) = EDetahat(1,1);
MOMvec(8) = EDetahat(2,1);
MOMvec(9) = EDetahat(3,1);
MOMvec(10) = EDetahat(1,2);
MOMvec(11) = EDetahat(2,2);
MOMvec(12) = EDetahat(3,2);
MOMvec(13) = EDetahat(1,3);
MOMvec(14) = EDetahat(2,3);
MOMvec(15) = EDetahat(3,3);
MOMvec(16) = EDetahat(1,4);
MOMvec(17) = EDetahat(2,4);
MOMvec(18) = EDetahat(3,4);

%do minimization of GMM objective
Aineq = zeros(3,Nparams);
Aineq(1,1) = -1;
Aineq(2,5)=-1;
Aineq(3,9)=-1;
Bineq = zeros(3,1);

boundval = 1.75; 

extraoutput = 0;
paramhat = fmincon(@fGMMobj,paramguess,Aineq,Bineq,[],[],-boundval*ones(Nparams,1),boundval*ones(Nparams,1),[],fminconopt);
minOBJ = fGMMobj(paramhat);
extraoutput=1;
momhat = fGMMobj(paramhat);
    
Bhat = zeros(NX,NX);
Bhat(1,1) = paramhat(1); 
Bhat(2,1) = paramhat(2);
Bhat(3,1) = paramhat(3);
Bhat(1,2) = paramhat(4); 
Bhat(2,2) = paramhat(5);
Bhat(3,2) = paramhat(6);
Bhat(1,3) = paramhat(7); 
Bhat(2,3) = paramhat(8);
Bhat(3,3) = paramhat(9);
Dcoeffhat=zeros(ND,2);
Dcoeffhat(1,1) = paramhat(10);
Dcoeffhat(2,1) = paramhat(11);
Dcoeffhat(3,1) = paramhat(12);
Dcoeffhat(4,1) = paramhat(13);
Dcoeffhat(1,2) = paramhat(14);
Dcoeffhat(2,2) = paramhat(15);
Dcoeffhat(3,2) = paramhat(16);
Dcoeffhat(4,2) = paramhat(17);


%implied shock decomposition
epshat = zeros(T,NX);
for t=1:T;
    etadum = etahat(t,:)';
    epsdum = Bhat\etadum;
    epshat(t,:) = epsdum';
end;


%IRFs
B1tilde = zeros(NX*Nlags,NX*Nlags);
B1tilde((1:NX),(1:(NX*Nlags)))=B1hat;
if (Nlags>1);
for ct=1:(Nlags-1);
    B1tilde(ct*NX+(1:NX),(ct-1)*NX+(1:NX))=eye(NX);
end;
end;
Btilde = zeros(NX*Nlags,NX);
Btilde(1:NX,1:NX)=Bhat;

IRF = zeros(lengthIRF,NX,NX);
for varct=1:NX;
for t=1:lengthIRF;
    IRFvec = (B1tilde^(t-1))*Btilde(:,varct);
    IRF(t,:,varct)= IRFvec(1:NX)';
end;


IRF(:,:,varct)=sqrt(var(X(:,varct)))*IRF(:,:,varct)/Bhat(varct,varct);


end;



IRF_S_TO_Y = IRF(1:15,1,3);
save IV_VAR.mat


