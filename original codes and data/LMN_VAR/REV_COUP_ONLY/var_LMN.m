%%%%%%%
% var_LMN.m
%
% This MATLAB code computes rev./coup. only admissible sets for the 
% disaster event restrictions VAR in 
% "Using Disasters to Estimate the Impact of Uncertainty"
% by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
%
% Initially run on MATLAB R2021a
%
% Questions to Stephen Terry
% stephenjamesterry@gmail.com
%%%%%%%

close all
clear all
clc

%set seed
rng(25041);

Nlags = 12; %number of VAR lags

N = 1500000; %number of draws of random matrices
Tirf = 50; %number of periods to plot
lwidnum = 2.5; %line width
fsizenum = 14; %font size

%read in estimation series & residuals from STATA
data = importdata('VARout.csv');
data = data.data;

%define events
pol_event = data(:,11); %this is a coup indicator
ter_event = data(:,12); %this is a terror attack indicator
nat_event = data(:,13); %this is a nat disaster indicator
rev_event = data(:,14); %this is revolution indicator

%read in and process coefficient vectors from STATA
Avec_y = importdata('Avec_y.txt');
Avec_y = Avec_y.data;
Avec_y = Avec_y(1:(3*Nlags));

Avec_first = importdata('Avec_first.txt');
Avec_first = Avec_first.data;
Avec_first = Avec_first(1:(3*Nlags));

Avec_second = importdata('Avec_second.txt');
Avec_second = Avec_second.data;
Avec_second = Avec_second(1:(3*Nlags));

%convert to raw coefficient matrices
Amats = zeros(3,3,Nlags);
for k=1:Nlags
    for varct=1:3
        Amats(1,varct,k) = Avec_y((varct-1)*Nlags+k);
        Amats(2,varct,k) = Avec_first((varct-1)*Nlags+k);
        Amats(3,varct,k) = Avec_second((varct-1)*Nlags+k);
    end
end

%put coeff mats in companion form
Abar = zeros(3*Nlags,3*Nlags);

for k=1:Nlags
   stct=(k-1)*3+1;
   endct = k*3;
   Abar(1:3,stct:endct) = squeeze(Amats(:,:,k));
end

for k=1:(Nlags-1)
   stct_x = 3 + (k-1)*3+1;
   endct_x = 3+k*3;
   stct_y = (k-1)*3+1;
   endct_y = k*3;
   Abar(stct_x:endct_x,stct_y:endct_y) = eye(3,3);
   
end


forecast_IRFbar = zeros(3*Nlags,3,Tirf);
forecast_Bbar = zeros(3*Nlags,3);
forecast_Bbar(1:3,1:3) = [1 0 0; 0 1 0; 0 0 1];
  for horz = 1:Tirf
     forecast_IRFbar(:,:,horz) = (Abar^(horz-1))*forecast_Bbar;
  end
  
IRFaxis = [ 1 6*4 -5 5];

%compute variance covariance matrix of the residuals
resids = data(:,6:8);
OMEGA = cov(resids);
SIGMA = chol(OMEGA,'lower');
stdev_shocks =sqrt(diag(OMEGA));

randB_store = zeros(3,3,N);
mean_shocks_store = zeros(3,N);
mean_shocks_store_coup = zeros(3,N);
mean_shocks_store_ter = zeros(3,N);
mean_shocks_store_nat = zeros(3,N);

for randct=1:N
    
    
    %generate random rotation
    randM = randn(3,3);
    [randQ,randR] = qr(randM);

    randB = SIGMA*randQ;
    for shockct=1:3
        if (randB(shockct,shockct)<0)
         randB(:,shockct) = -randB(:,shockct);
        end
    end
    
    randB_store(:,:,randct) = randB;
      
    %compute implied shocks
    shocks = randB\resids';
    shocks = shocks';
    
    %mean shocks
    mean_shocks = mean(shocks(rev_event==1,:),1);
    mean_shocks_store(:,randct) = mean_shocks;
    
    mean_shocks_coup = mean(shocks(pol_event==1,:),1);
    mean_shocks_store_coup(:,randct) = mean_shocks_coup;
    
    mean_shocks_ter = mean(shocks(ter_event==1,:),1);
    mean_shocks_store_ter(:,randct) = mean_shocks_ter;
    
    mean_shocks_nat = mean(shocks(nat_event==1,:),1);
    mean_shocks_store_nat(:,randct) = mean_shocks_nat;
        
        
    if (mod(randct,100)==1)
       disp(['Done with ' num2str(100*randct/N) '%.']) 
    end
    
    
end

%now, process the candidates to see which IRFs are admissable
B_admit_lb = zeros(3,3)+1000000;
B_admit_ub = zeros(3,3)-1000000;
IRF_admit_lb = zeros(3,3,Tirf)+1000000;
IRF_admit_ub = zeros(3,3,Tirf)-1000000;
IRF_store = zeros(3,3,Tirf,N);
IRF_maxg = zeros(3,3,Tirf);
admit_vec = zeros(N,1);
admit_ind = [];

Nadmit= 0 ;
for randct=1:N
      

                 if (mean_shocks_store(3,randct)>0.15)&...
           (mean_shocks_store(2,randct)<-0.1)&...
           (mean_shocks_store_coup(3,randct)>0.15)&...
           (mean_shocks_store_coup(2,randct)>0.1);

     admit_vec(randct) = 1;
     admit_ind = [admit_ind; randct];
    Nadmit = Nadmit +1;

    %compute implied IRF
    rand_Bbar = zeros(3*Nlags,3);
    rand_Bbar(1:3,1:3) = randB_store(:,:,randct);
    rand_IRFbar = zeros(3*Nlags,3,Tirf);
    for horz = 1:Tirf
       rand_IRFbar(:,:,horz) = (Abar^(horz-1))*rand_Bbar;
    end

        
        
    rand_IRF = rand_IRFbar(1:3,1:3,1:Tirf);


    IRF_store(:,:,:,randct) = rand_IRF;



    for varct=1:3
       for shockct=1:3
             B_admit_lb(varct,shockct) = min([B_admit_lb(varct,shockct) randB_store(varct,shockct,randct)]);
             B_admit_ub(varct,shockct) = max([B_admit_ub(varct,shockct) randB_store(varct,shockct,randct)]);
             for t=1:Tirf
                IRF_admit_lb(varct,shockct,t) = min([IRF_admit_lb(varct,shockct,t) rand_IRF(varct,shockct,t)]);
                IRF_admit_ub(varct,shockct,t) = max([IRF_admit_ub(varct,shockct,t) rand_IRF(varct,shockct,t)]);
             end
       end
    end
    
    end
    
end

IRF_med = zeros(3,3,Tirf);
for varct=1:3
    for shockct=1:3
        for t=1:Tirf
            IRF_med(varct,shockct,t) = median(squeeze(IRF_store(varct,shockct,t,admit_ind)));
        end
    end
end



IMPACT_HIST = squeeze(IRF_store(1,3,1,admit_ind));


clear IRF_store

save LMN_WORK.mat;

cd ..