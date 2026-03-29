program VOL_GROWTH_wrapper
use base_lib
use omp_lib
implicit none

!!!!!!
!VOL_GROWTH_wrapper.f90
!
!This code performs the optimization of the indirect inference objective
!function detailed in "Using Disasters to Estimate the Impact of Uncertainty,"
! by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry
!
!Questions to Stephen Terry
!stephenjamesterry@gmail.com
!!!!!!

!!!!DECLARATION BLOCK
integer, parameter :: nummom = 20
integer, parameter :: numparam = 8
integer, parameter :: maxevals = 1000

integer :: evalct,paramct,seedintout,newrun,npart,xquicknum,maxpsoit,seeddimout,ct,psoseed
integer, allocatable :: seedarrayout(:)

double precision :: x(numparam),GMMout,lb(numparam),ub(numparam),randparamshocks(maxevals,numparam),&
	randval,xpsotol,fpsotol,xquicktol,phipso(2),xpso(numparam),fpso

!PSO setup and newrun flag
newrun = 1
psoseed=8791
npart = 75
xpsotol = 1.0e-3
fpsotol = 1.0e-3
xquicktol = 1.0e-3
xquicknum = 5
phipso = (/2.05,2.05/)
maxpsoit = 5000

!indexing for parameter vector x
!1 - levels impact of nat dis
!2 - levels impact of pol shock
!3 - levels impact of revolution
!4 - levels impact of terrorist
!5 - unc impact of nat dis
!6 - unc impact of pol shock
!7 - unc impact of revolution
!8 - unc impact of terrorist

!best fit overall from initial uniform exploration
x(1) =  -0.034
x(2) =  0.054
x(3) =  -41.486
x(4) =  -3.950
x(5) =  0.014
x(6) =  0.853
x(7) =  0.816
x(8) =  0.110

GMMout =fGMM(x)

!set upper and lower boundaries for exploration
lb(1) = -1.5; ub(1) = 0.1 ! nat dis levels
lb(2) = -0.1; ub(2) = 0.25 ! pol shock levels
lb(3) = -75.0; ub(3) = -35.5 ! revolution levels
lb(4) = -4.5; ub(4) = -2.5 ! terrorist levels

lb(5) = 0.001; ub(5) = 0.35 ! nat dis unc, in prob units
lb(6) = 0.5; ub(6) = 0.999 ! pol shock unc, in prob units
lb(7) = 0.75; ub(7) = 0.999 ! revolution unc, in prob units
lb(8) = 0.001; ub(8) = 0.25 ! terrorist unc, in prob units


open(8,file="lb_ub.txt")
do paramct=1,numparam
write(8,*) lb(paramct),ub(paramct)
end do !paramct
close(8)

!do the PSO optimization
call psorestart(x,GMMout,fGMM,lb,ub,numparam,npart,xpsotol,xquicktol,xquicknum,fpsotol,maxpsoit,&
    phipso,newrun,psoseed)

!!!!INTERNAL FUNCTION AND SUBROUTINES BLOCK
contains

double precision function fGMM(x)
implicit none

double precision :: x(8) !input vector with the parameters controlling the model first stage

integer :: znum,anum,snum,knum,lnum,numstates,numexog,numendog,numper,numdiscard,seedint,vfmaxit,&
    accelmaxit,numconstants,ct,zct,act,sct,kct,lmin1ct,zprimect,aprimect,sprimect,primect,endogct,&
    exogct,polct,vfct,accelct,polstar,exogprimect,t,numsimIRF,lengthIRF,shockperIRF,numperIRF,&
    numdiscIRF,ainit,sinit,kbarnum,numfcst,fcstct,lct,kprimect,ind,piter,maxpit,doVFI,distinit,&
    GEct,maxGEit,perct,doIRF,pcutoff,checkbounds,shocklengthIRF,simct,GEerrorswitch,exitflag,&
    singleshock,smin1ct,seeddim,pnum,pind,nfirms,nfirmspub,firmct,zinit,numyr,Ncountries,Tper,&
    countryct,waitsec,maxwaits

integer, allocatable :: exog_pos(:,:),endog_pos(:,:),loopind(:,:),polmat(:,:,:),polmatold(:,:,:),&
    asimpos(:),ssimpos(:),asimposIRF(:,:),ssimposIRF(:,:),kbarfcstinds(:,:,:,:),polmatp(:,:),&
    kprime_posp(:,:),lpol_posp(:,:),piterstoresim(:),piterstoresimIRF(:,:),seedarray(:),&
    polmatp_interp(:,:,:),endogfirmpos(:,:),zfirmpos(:,:),firmpos(:,:),&
	DISASTERpos(:,:)

double precision :: start,finish,alpha,nu,theta,deltak,deltan,beta,capirrev,capfix,hirelin,firelin,&
    labfix,pval,wval,ajump,zjump,uncpers,uncfreq,rhoz,sigmaz,rhoa,sigmaa,zmin,zmax,amin,amax,kmin,kmax,&
    lmin,lmax,vferrortol,zval,aval,sval,kval,lmin1val,kprimeval,lval,ival,yval,ackval,aclval,vferror,&
    polerror,hval,changetol,kbarmin,kbarmax,kbarval,kbarfcstval,pfcstval,wfcstval,Vnextval,weight,&
    perror,perrortol,Yvalp,Ivalp,ACkvalp,AClvalp,Cvalp,pupdate,Kbarprimevalp,Hvalp,Lvalp,GEupdate,&
    xmean,ymean,x2mean,xymean,xval,kfcsterror,pfcsterror,fcsterrortol,pvala,pvalb,pvalc,fa,fb,fc,plb,&
    pub,pwidth,disttol,pwindow,KRMSEchange,pRMSEchange,KR2change,pR2change,RMSEchangetol,R2changetol,&
    maxDenHaanchangetol,avgDenHaanchangetol,pmaxDenHaanchange,KmaxDenHaanchange,pavgDenHaanchange,&
    KavgDenHaanchange,pwgt,returnmean,returnstdev,DISASTERprobs(4),DISASTERuncprobs(4),meanshifta,&
    DISASTERlev(4),meanval,sdval,FIRST_STAGE_MACRO_DATA(4,2),FIRST_STAGE_MICRO_DATA(4,2),&
    SECOND_STAGE_MACRO_DATA(2),SECOND_STAGE_MICRO_DATA(2),MATLAB_MOMS(20),&
    DATA_MOMS(20),GMMobj,totweightval,weighttol,firstsecondprob,&
		highuncerg,DATA_SE(20),aggretmean,aggretmean2,aggretstd,aggretMEmult

double precision, allocatable :: V(:,:,:),asimshocks(:),ssimshocks(:),pr_mat_z(:,:,:),&
    pr_mat_a(:,:,:),pr_mat_s(:,:),pr_mat(:,:),z0(:),k0(:),l0(:),a0(:),sigmazgrid(:),sigmaagrid(:),&
    MATLABconstants(:),exog_key(:,:),endog_key(:,:),EVmat(:,:,:),Vold(:,:,:),RHSvec(:),kprime_key(:,:,:),&
    lpol_key(:,:,:),kprime_pos(:,:,:),lpol_pos(:,:,:),distzkl(:,:,:),Ysim(:),Ksim(:),Lsim(:),Isim(:),Hsim(:),&
    ACksim(:),AClsim(:),TOPksim(:),BOTksim(:),TOPlsim(:),BOTlsim(:),ADJksim(:),ADJlsim(:),asimshocksIRF(:,:),&
    ssimshocksIRF(:,:),YsimIRF(:,:),KsimIRF(:,:),LsimIRF(:,:),IsimIRF(:,:),HsimIRF(:,:),ACksimIRF(:,:),&
    AClsimIRF(:,:),ADJksimIRF(:,:),ADJlsimIRF(:,:),kbar0(:),kbarfcstmat(:,:,:,:),pfcstmat(:,:,:,:),kfcstcoeff(:,:,:,:),&
    pfcstcoeff(:,:,:,:),wfcstmat(:,:,:,:),Ymat(:,:,:,:),Imat(:,:),ACkmat(:,:,:,:,:),AClmat(:,:,:,:,:),&
    WLmat(:,:,:,:,:),kbarfcstweights(:,:,:,:),psim(:),EVmatp(:,:),Csim(:),Kfcstsim(:),pfcstsim(:),&
    kcoeffstore(:,:,:,:,:),pcoeffstore(:,:,:,:,:),kfcstcoeffnew(:,:,:,:),pfcstcoeffnew(:,:,:,:),&
    kfcsterrorstore(:),pfcsterrorstore(:),psimIRF(:,:),CsimIRF(:,:),KfcstsimIRF(:,:),pfcstsimIRF(:,:),&
    kfcsterrorstoreIRF(:,:),pfcsterrorstoreIRF(:,:),distzklbefore(:,:),distzklafter(:,:),pR2store(:,:,:,:),&
    KR2store(:,:,:,:),pRMSEstore(:,:,:,:),KRMSEstore(:,:,:,:),pSSE(:,:,:),KSSE(:,:,:),pSST(:,:,:),KSST(:,:,:),Kmean(:,:,:),&
    pmean(:,:,:),KRMSEchangestore(:),pRMSEchangestore(:),KR2changestore(:),pR2changestore(:),pDenHaanfcst(:),&
    KDenHaanfcst(:),pmaxDenHaanstat(:),KmaxDenHaanstat(:),pavgDenHaanstat(:),KavgDenHaanstat(:),&
    pmaxDenHaanchangestore(:),KmaxDenHaanchangestore(:),pavgDenHaanchangestore(:),KavgDenHaanchangestore(:),p0(:),&
    ep0(:),Cp0(:),Yp0(:),Ip0(:),ACkp0(:),AClp0(:),Kbarprimep0(:),Hp0(:),Lp0(:),vmatp_interp(:,:,:),&
    yfirm(:,:),vfirm(:,:),zfirmshocks(:,:),pfirmshocks(:,:),FIRSTsim(:),SECONDsim(:),GDPsim(:),GROWTHsim(:),&
	normfirmnoise(:,:),returnfirm(:,:),returnfirmnoise(:,:),GROWTHsimYR(:),FIRSTsimYR(:),SECONDsimYR(:),&
	DISASTERshocks(:,:),achgshocks(:),lfirm(:,:),schgshocks(:),returnfirmsd(:,:),dfirm(:,:),returnfirmsdnoise(:,:),&
	ascorrshocks(:),ushocks1(:,:),ushocks2(:,:),aggretnoise(:),normaggretnoise(:,:),ushocksagg1(:,:),ushocksagg2(:,:)

character*200 :: matlabstr !this is the string variable used to call MATLAB below


character*200 :: outputfilestr !this is the string variable used to call MATLAB below

!!!!END DECLARATION BLOCK

!these two vectors control size of first stage within the model

DISASTERlev(1) = x(1)
DISASTERlev(2) = x(2)
DISASTERlev(3) = x(3)
DISASTERlev(4) = x(4)
DISASTERuncprobs(1) = x(5)
DISASTERuncprobs(2) = x(6)
DISASTERuncprobs(3) = x(7)
DISASTERuncprobs(4) = x(8)

aggretMEmult = 0.0
firstsecondprob = 0.0 !prob of moving to lowest prod when unc is high


write(*,*) "###########################################################"
write(*,*) "###########################################################"
write(*,*) "Starting IV Matching Function Evaluation"

write(*,*) "Disaster Levels Params"
do ct=1,4
write(*,*) DISASTERlev(ct)
end do !ct
write(*,*) "Disaster Unc Params"
do ct=1,4
write(*,*) DISASTERuncprobs(ct)
end do !ct


!!!!PROGRAM SETUP BLOCK

doVFI = 1 !1 implies do VFI, 0 implies read from file
distinit = 0; !1 implies read initial distribution from file (only matters for GE solution, not IRF)
doIRF = 0; !1 implies do the IRF after GE iteration is complete
singleshock = 1; !1 implies a single shock period for the IRF, 0 implies a low unc period followed by high unc
checkbounds = 1; !1 implies that after the uncond sim you will check to see how often you the top of the grid

!GE convergence criteria
!1 implies convergence of coefficients
!2 implies convergence of RMSE in changes
!3 implies convergence of R^2 in changes
!4 implies convergence of max Den Haan stat in changes
!5 implies convergence of avg Den Haan stat in changes
GEerrorswitch = 4;

!call the clock to get the start time and start writing to the log file
start = omp_get_wtime()

!!!!END PROGRAM SETUP BLOCK

!!!!!PARAMETER ENTRY BLOCK

!indexing setup for MATLAB_MOMS vector
! 1 - ann, nat dis first stage for first mom
! 2 - ann, pol shock first stage for first mom
! 3 - ann, revolution first stage for first mom
! 4 - ann, terror first stage for first mom
! 5 - ann, nat dis first stage for second mom
! 6 - ann, pol shock first stage for second mom
! 7 - ann, revolution first stage for second mom
! 8 - ann, terror first stage for second mom
! 9 - ann, first moment coeff OLS
! 10 - ann, second moment coeff OLS
! 11 - ann, first moment coeff second stage
! 12 - ann, second moment coeff second stage

!enter first & second stages from the data (MACRO SAMPLE)
FIRST_STAGE_MACRO_DATA(1,:) = (/ -0.071, -0.028 /) !Nat Dis
FIRST_STAGE_MACRO_DATA(2,:) = (/ 1.657, 1.693 /) !Pol Shock
FIRST_STAGE_MACRO_DATA(3,:) = (/ -6.154, 7.841 /) !Revolution
FIRST_STAGE_MACRO_DATA(4,:) = (/ -0.047, -0.011 /) !Terrorist

SECOND_STAGE_MACRO_DATA = (/1.557, -3.859/) !First Moment, Second Moment

FIRST_STAGE_MICRO_DATA(1,:) = (/ -0.147, 0.004/) !Nat Dis
FIRST_STAGE_MICRO_DATA(2,:) = (/ 1.852, 0.508 /) !Pol Shock
FIRST_STAGE_MICRO_DATA(3,:) = (/ -4.818, 3.201 /) !Revolution
FIRST_STAGE_MICRO_DATA(4,:) = (/ -0.117, 0.133 /) !Terrorist

SECOND_STAGE_MICRO_DATA = (/0.736, -9.735/) !First Moment, Second Moment


!insert data moments into combined vector

!data moment organization
!1 - nat dis to first, macro
!2 - pol shock to first, macro
!3 - revolution to first, macro
!4 - terror to first, macro
!5 - nat dis to second, macro
!6 - pol shock to second, macro
!7 - revolution to second, macro
!8 - terror to second, macro
!9 - second stage, first moment, macro
!10 - second stage, second moment, macro
!11 - nat dis to first, micro
!12 - pol shock to first, micro
!13 - revolution to first, micro
!14 - terror to first, micro
!15 - nat dis to second, micro
!16 - pol shock to second, micro
!17 - revolution to second, micro
!18 - terror to second, micro
!19 - second stage, first moment, micro
!20 - second stage, second moment, micro

DATA_MOMS(:) = 0.0
DATA_MOMS(1:4) = FIRST_STAGE_MACRO_DATA(:,1)
DATA_MOMS(5:8) = FIRST_STAGE_MACRO_DATA(:,2)
DATA_MOMS(9:10) = SECOND_STAGE_MACRO_DATA
DATA_MOMS(11:14) = FIRST_STAGE_MICRO_DATA(:,1)
DATA_MOMS(15:18) = FIRST_STAGE_MICRO_DATA(:,2)
DATA_MOMS(19:20) = SECOND_STAGE_MICRO_DATA

open(8,file="DATA_MOMS.txt")
do ct=1,20
write(8,*) DATA_MOMS(ct)
end do !ct
close(8)

!insert the precision of each moment
DATA_SE(:) = 0.0

!first stage, levels LHS, macro
DATA_SE(1) = 0.1059
DATA_SE(2) = 0.0551
DATA_SE(3) = 1.0835
DATA_SE(4) = 0.0514

!first stage, vol LHS, macro
DATA_SE(5) = 0.0823
DATA_SE(6) = 0.1159
DATA_SE(7) = 2.2358
DATA_SE(8) = 0.0490

!second stage, growth LHS, macro
DATA_SE(9) = 0.2906
DATA_SE(10) = 0.2844

!first stage, levels LHS, micro
DATA_SE(11) = 0.112
DATA_SE(12) = 0.085
DATA_SE(13) = 1.198
DATA_SE(14) = 0.044

!first stage, vol LHS, micro
DATA_SE(15) = 0.102
DATA_SE(16) = 0.130
DATA_SE(17) = 1.275
DATA_SE(18) = 0.083

!second stage, growth LHS, micro
DATA_SE(19) = 0.558
DATA_SE(20) = 1.533

!(Nat Dis, Pol Shock, Revolution, Terrorist)
!data frequencies
DISASTERprobs = (/0.242,0.03,0.011,0.008/)

!technology and adjustment costs
alpha = 0.25
nu = 0.5
theta = 2.0
deltak = 0.026
deltan = 0.088
beta = 0.95 ** 0.25
capirrev = 0.339;
capfix = 0.0
hirelin = 0.018*4.0;
firelin = hirelin
labfix = 0.024*4;

!grid sizes for exog processes
znum = 9; !idio prod
anum = 21; !agg prod
snum = 2; !volatility

!grid size for interpolation of excess demand function
pnum = 1
pval = 1.34
plb = pval;
pub = pval

!grid sizes for endog processes
knum = 91
lnum = 37

knum = 150
lnum = 75

!grid size for the aggregate moments
kbarmin = 3.0
kbarmax = 10.0
kbarnum = 2

!set up the unc process
ajump = 1.60569110682638; !multiple of vol in agg case
zjump = 4.11699578856773; !multiple of vol in idio case

uncpers = 0.940556523390567; !cond. unc shock persistence
uncfreq = 0.0257892725462263; !cond. prob of unc shock

!set up the idio prod process
rhoz = 0.95;
sigmaz = 0.0507515557155377
zmin = exp( - 2.5 * ( ( sigmaz ** 2.0 ) / ( 1.0 - rhoz ** 2.0 ) ) ** 0.5 );
zmax = exp( 2.5 * ( ( sigmaz ** 2.0 ) / ( 1.0 - rhoz ** 2.0 ) ) ** 0.5);

!set up the agg prod process
rhoa = 0.95;
sigmaa = 0.00668420914017636
amin = 0.8
amax = 1.2

!determine the upward shift in the mean agg prod to undo agg disaster mean impact
highuncerg = uncfreq / (1.0+uncfreq-uncpers)
meanshifta = 0.0
do ct=1,4
	meanshifta = meanshifta + sigmaa*DISASTERprobs(ct)*DISASTERlev(ct)
end do !ct
meanshifta = meanshifta + highuncerg*firstsecondprob*log(amin)
meanshifta = dble(-1.0)*meanshifta


!adjust uncfreq down to compensate for frequency of unc shocks
do ct=1,4
	uncfreq = uncfreq - DISASTERprobs(ct)*DISASTERuncprobs(ct)
end do !ct


!set up the idio capital grid, respecting depreciation
kmin = 0.75
kmax = exp(log(kmin) - dble(knum-1) * log(1.0-deltak))

!set up the idio labor grid, respecting depreciation
lmin = 0.02
lmax = exp(log(lmin) - dble(lnum-1) * log(1.0-deltan))


!figure out total number of states
numexog = znum*anum*snum*snum
numendog = knum*lnum
numfcst = kbarnum
numstates = numexog * numendog * numfcst

!control the uncondtional simulation
Ncountries = 500
Tper = 100
numdiscard = 500
numper = Ncountries*Tper+numdiscard
seedint = 2501; call random_seed(size=seeddim)
ainit = 3; !initial agg prod val (uncond sim and IRF)
sinit = 1; !initial unc val (uncond sim only, 1 = low, 2 = high)

!control number of firms and public firms
nfirms = 800 !number of firms in GDP data
nfirmspub = 200 !number of firms in stock market data

zinit = 3 !initial firm productivity position

!control the IRF
numsimIRF = 2500; !number of shocked economies to simulate
lengthIRF = 100; !length of each economy
shockperIRF = 45; !period in which to shock the economy
shocklengthIRF = 5; !number of periods to have low unc
numdiscIRF = 45; !number of periods to discard in the IRF

numperIRF = numsimIRF * lengthIRF; !total number of periods in IRF simulation

!control the vf process
vfmaxit=50; !number of VF iterations
vferrortol = 1e-4; !VF error tolerance
accelmaxit = 200; !number of Howard acceleration iterations

!control the price-clearing process
maxpit = 50
perrortol = 1.0e-3
!plb = 1.0; !price lower boundary (initial bisection lb)
!pub = 1.8; !price upper boundary (initial bisection ub)
disttol = 1e-4; !tolerance for ignoring this point in the dist
pwindow = 0.025
pcutoff = 15

!tolerance for sum total of weight on the grid endpoints over the simulation at which to discard with GMMobj = 1000000
weighttol = 1e-4

!control the GE fcst rule update process
maxGEit = 1
GEupdate = 0.33; !step speed in fcst rule fixed-point iteration process
fcsterrortol = 0.0125; !the convergence tolerance for the forecast rules
RMSEchangetol = 0.001; !the convergence tolerance for RMSE changes
R2changetol = 0.01; !the convergence tolerance for R2 changes
maxDenHaanchangetol = 0.01; !the convergence tolerance for Den Haan max statistic
avgDenHaanchangetol = 0.01; !the convergence tolerance for Den Haan avg statistic

!control what is considered a change
changetol = 1e-10; !what is a change for the purposes of AC functions?

!how many parameters to pass to MATLAB
numconstants = 20; !this is relevant for input/output only in pass to MATLAB

waitsec = 2 !the number of seconds between checks for MATLAB completion
maxwaits = 20 !the max number of waits to check if the file exists


!!!!END PARAMETER ENTRY BLOCK

!!!!SOME INITIALIZATIONS/GRID SETUPS

!do allocations for all variable-dimension arrays
allocate(k0(knum),l0(lnum),z0(znum),a0(anum),pr_mat_z(znum,znum,snum),pr_mat_a(anum,anum,snum),pr_mat_s(snum,snum),&
    sigmaagrid(snum),sigmazgrid(snum),MATLABconstants(numconstants),exog_key(numexog,4),exog_pos(numexog,4),endog_key(numendog,2),&
    endog_pos(numendog,2),V(numendog,numexog,numfcst),pr_mat(numexog,numexog),EVmat(numexog,numendog,numfcst),&
    loopind(numstates,3),polmat(numendog,numexog,numfcst),Vold(numendog,numexog,numfcst),polmatold(numendog,numexog,numfcst),&
    RHSvec(numendog),kprime_key(numendog,numexog,numfcst),lpol_key(numendog,numexog,numfcst),&
    kprime_pos(numendog,numexog,numfcst),lpol_pos(numendog,numexog,numfcst),asimpos(numper),ssimpos(numper),&
    distzkl(znum,numendog,numper),&
    asimshocks(numper),ssimshocks(numper),Ysim(numper),Ksim(numper),Lsim(numper),Isim(numper),Hsim(numper),ACksim(numper),&
    AClsim(numper),TOPksim(numper),BOTksim(numper),TOPlsim(numper),BOTlsim(numper),ADJksim(numper),ADJlsim(numper),psim(numper),&
    asimposIRF(lengthIRF,numsimIRF),ssimposIRF(lengthIRF,numsimIRF),asimshocksIRF(lengthIRF,numsimIRF),&
    ssimshocksIRF(lengthIRF,numsimIRF),&
    YsimIRF(lengthIRF,numsimIRF),KsimIRF(lengthIRF,numsimIRF),LsimIRF(lengthIRF,numsimIRF),IsimIRF(lengthIRF,numsimIRF),&
    HsimIRF(lengthIRF,numsimIRF),ACksimIRF(lengthIRF,numsimIRF),AClsimIRF(lengthIRF,numsimIRF),ADJksimIRF(lengthIRF,numsimIRF),&
    ADJlsimIRF(lengthIRF,numsimIRF),kbar0(kbarnum),&
    kbarfcstmat(anum,snum,snum,numfcst),pfcstmat(anum,snum,snum,numfcst),kfcstcoeff(anum,snum,snum,2),pfcstcoeff(anum,snum,snum,2),&
    wfcstmat(anum,snum,snum,numfcst),Ymat(znum,anum,knum,lnum),Imat(knum,knum),ACkmat(znum,anum,knum,lnum,knum),&
    AClmat(numexog,knum,lnum,lnum,numfcst),WLmat(anum,snum,snum,numfcst,lnum),kbarfcstweights(anum,snum,snum,numfcst),&
    kbarfcstinds(anum,snum,snum,numfcst),EVmatp(znum,numendog),polmatp(znum,numendog),kprime_posp(znum,numendog),&
    lpol_posp(znum,numendog),Csim(numper),Kfcstsim(numper),pfcstsim(numper),kcoeffstore(anum,snum,snum,2,maxGEit),&
    pcoeffstore(anum,snum,snum,2,maxGEit),kfcstcoeffnew(anum,snum,snum,2),pfcstcoeffnew(anum,snum,snum,2),piterstoresim(numper),&
    kfcsterrorstore(maxGEit),pfcsterrorstore(maxGEit),psimIRF(lengthIRF,numsimIRF),CsimIRF(lengthIRF,numsimIRF),&
    KfcstsimIRF(lengthIRF,numsimIRF),&
    pfcstsimIRF(lengthIRF,numsimIRF),kfcsterrorstoreIRF(lengthIRF,numsimIRF),pfcsterrorstoreIRF(lengthIRF,numsimIRF),&
    piterstoresimIRF(lengthIRF,numsimIRF),&
    distzklbefore(znum,numendog),distzklafter(znum,numendog),pR2store(anum,snum,snum,maxGEit),KR2store(anum,snum,snum,maxGEit),&
    pRMSEstore(anum,snum,snum,maxGEit),KRMSEstore(anum,snum,snum,maxGEit),&
    pSSE(anum,snum,snum),KSSE(anum,snum,snum),pSST(anum,snum,snum),&
    KSST(anum,snum,snum),Kmean(anum,snum,snum),pmean(anum,snum,snum),KRMSEchangestore(maxGEit),pRMSEchangestore(maxGEit),&
    KR2changestore(maxGEit),pR2changestore(maxGEit),pDenHaanfcst(numper),KDenHaanfcst(numper),pmaxDenHaanstat(maxGEit),&
    KmaxDenHaanstat(maxGEit),pavgDenHaanstat(maxGEit),KavgDenHaanstat(maxGEit),pmaxDenHaanchangestore(maxGEit),&
    KmaxDenHaanchangestore(maxGEit),pavgDenHaanchangestore(maxGEit),KavgDenHaanchangestore(maxGEit),&
    seedarray(seeddim),p0(pnum),ep0(pnum),Cp0(pnum),Yp0(pnum),Ip0(pnum),ACkp0(pnum),AClp0(pnum),&
    Kbarprimep0(pnum),Hp0(pnum),Lp0(pnum),polmatp_interp(znum,numendog,pnum),vmatp_interp(znum,numendog,pnum),&
    endogfirmpos(numper,nfirms),zfirmpos(numper,nfirms),yfirm(numper,nfirms),vfirm(numper,nfirms),zfirmshocks(numper,nfirms),&
    pfirmshocks(numper,nfirms),FIRSTsim(numper),SECONDsim(numper),GDPsim(numper),GROWTHsim(numper),firmpos(znum,numendog),&
    normfirmnoise(numper,nfirms),ushocks1(numper,nfirms),ushocks2(numper,nfirms),returnfirm(numper,nfirms),&
		returnfirmnoise(numper,nfirms),GROWTHsimYR(numper),FIRSTsimYR(numper),SECONDsimYR(numper),DISASTERpos(numper,4),&
		DISASTERshocks(numper,4),achgshocks(numper),lfirm(numper,nfirms),schgshocks(numper),returnfirmsd(numper,nfirms),&
		dfirm(numper,nfirms),returnfirmsdnoise(numper,nfirms),ascorrshocks(numper),aggretnoise(numper),&
		normaggretnoise(numper,1),ushocksagg1(numper,1),ushocksagg2(numper,1))

!set up the idiosyncratic capital and labor grids, using the subroutine linspace which replicates MATLAB function
call linspace(k0,log(kmin),log(kmax),knum); k0=exp(k0);
call linspace(l0,log(lmin),log(lmax),lnum); l0=exp(l0);

p0 = (/pval/)
kbar0 = (/kbarmin,kbarmax/)

!set up the sigma or volatility grids and transition matrix
sigmazgrid = (/ sigmaz , zjump*sigmaz /)
sigmaagrid = (/ sigmaa , ajump*sigmaa /)
pr_mat_s(1,:) = (/1.0 - uncfreq , uncfreq /); pr_mat_s(1,:) = pr_mat_s(1,:)/sum(pr_mat_s(1,:))
pr_mat_s(2,:) = (/1.0 - uncpers , uncpers /); pr_mat_s(2,:) = pr_mat_s(2,:)/sum(pr_mat_s(2,:))

kbarfcstmat(:,:,:,:) = (kbar0(1)+kbar0(2))/dble(2.0)
pfcstmat(:,:,:,:) = pval
wfcstmat(:,:,:,:) = theta/pval
kbarfcstinds(:,:,:,:) = 1
kbarfcstweights(:,:,:,:) = 0.0

!compute the idio & aggr transition matrices, via Tauchen (1986) adapted to stochastical vol case
call unctauchen(pr_mat_z,z0,znum,zmin,zmax,snum,sigmazgrid,rhoz,dble(0.0))
call unctauchen(pr_mat_a,a0,anum,amin,amax,snum,sigmaagrid,rhoa,meanshifta)

!index the exogenous variables (z,a,s) & create the unified transition matrix over exogenous variables
ct=0
do zct=1,znum; do act=1,anum; do sct=1,snum; do smin1ct=1,snum;
    ct=ct+1;

    !indexing the variables
    exog_pos(ct,1) = zct; exog_key(ct,1) = z0(zct); !insert z
    exog_pos(ct,2) = act; exog_key(ct,2) = a0(act); !insert a
    exog_pos(ct,3) = sct; exog_key(ct,3) = dble(sct); !insert s
    exog_pos(ct,4) = smin1ct; exog_key(ct,4) = dble(smin1ct); !insert s_{-1}

    !creating the transition matrix
    do zprimect=1,znum; do aprimect=1,anum; do sprimect=1,snum
        !note that this only indexes numexog/snum of the next period's values, since sct today is smin1ct tomorrow and is FIXED!
        primect=(zprimect-1)*anum*snum*snum + (aprimect-1)*snum*snum + (sprimect-1)*snum + sct
        pr_mat(ct,primect) = pr_mat_z(zct,zprimect,sct) * pr_mat_a(act,aprimect,sct) * pr_mat_s(sct,sprimect)
    end do; end do; end do;

    !correction for roundoff error
    pr_mat(ct,:) = pr_mat(ct,:) / sum( pr_mat(ct,:) )

end do; end do; end do; end do;

!call the random numbers for unconditional simulation - occurs outside of fcst loop
!this werid block seeds the draws
do ct=1,seeddim
    seedarray(ct) = seedint+ct
end do !ct
call random_seed(put=seedarray)
call random_number(asimshocks); !U(0,1) in each cell of asimshocks
call random_number(ssimshocks); !U(0,1) in each cell of ssimshocks
call random_number(DISASTERshocks) !U(0,1) in each cell of DISASTERshocks
call random_number(achgshocks) !U(0,1) in each cell of achgshocks
call random_number(schgshocks) !U(0,1) in each cell of schgshocks
call random_number(ascorrshocks) !U(0,1) in each cell of ascorrshocks

!based on the random shocks drawn, simulate the aggregate exogenous processes, occurs outside of fcst loop
call uncexogsim(ssimshocks,asimshocks,DISASTERshocks,DISASTERprobs,&
	ssimpos,asimpos,pr_mat_s,pr_mat_a,ainit,sinit,snum,anum,numper,a0,achgshocks,DISASTERpos,sigmaa,&
	schgshocks,DISASTERuncprobs,DISASTERlev,ascorrshocks,firstsecondprob)

!draw random shocks for the firm simulation
call random_number(zfirmshocks)
call random_number(pfirmshocks)

!then, based on the zfirmshocks matrix as well as ssimpos, simulate individual firm productivity processes
call firmexogsim(zfirmshocks,zfirmpos,ssimpos,numper,nfirms,pr_mat_z,znum,zinit,snum)

!draw uniform shocks for firm noise, and convert these to normal shocks
call random_number(ushocks1)
call random_number(ushocks2)
call boxmuller(normfirmnoise,numper,nfirms,ushocks1,ushocks2)


call random_number(ushocksagg1)
call random_number(ushocksagg2)
call boxmuller(normaggretnoise,numper,1,ushocksagg1,ushocksagg2)

!index the endog variables (k,l_{-1}), noting that the index is used for states and policies
ct=0
do kct=1,knum; do lmin1ct=1,lnum
    ct=ct+1
    endog_pos(ct,1) = kct; endog_key(ct,1) = k0(kct); !insert k
    endog_pos(ct,2) = lmin1ct; endog_key(ct,2) = l0(lmin1ct); !insert l_{-1}
end do; end do;

!indexing for parallelization over V: parallel loops will go over "ct," with
!each thread pulling out endogct, exogct, and fcstct
ct=0
do endogct=1,numendog
do exogct=1,numexog
do fcstct=1,numfcst
    ct=ct+1
    loopind(ct,1) = endogct; loopind(ct,2) = exogct; loopind(ct,3) = fcstct
end do !numfcst
end do !exogct
end do !endogct

!determine and write constants
MATLABconstants = (/ dble(znum),dble(knum),dble(lnum),dble(anum),dble(snum),&
    dble(kbarnum),dble(numper),dble(numdiscard),dble(numendog),dble(numexog),dble(numfcst),&
    dble(numstates),deltak,deltan,dble(numsimIRF),dble(lengthIRF),dble(shockperIRF),dble(numdiscIRF),&
    dble(shocklengthIRF),dble(singleshock)/)

open(8,file="MATLABconstants.txt")
do ct=1,numconstants
write(8,*) MATLABconstants(ct)
end do !ct
close(8)

!write labor grid
open(8,file="l0.txt")
do ct=1,lnum
write(8,*) l0(ct)
end do
close(8)

!write idio capital grid
open(8,file="k0.txt")
do ct=1,knum
write(8,*) k0(ct)
end do
close(8)

!write aggregate capital grid
open(8,file="kbar0.txt")
do ct=1,kbarnum
write(8,*) kbar0(ct)
end do
close(8)

!write idio prod grid
open(8,file="z0.txt")
do ct=1,znum
write(8,*) z0(ct)
end do
close(8)

!write agg prod grid
open(8,file="a0.txt")
do ct=1,anum
write(8,*) a0(ct)
end do
close(8)

!write endogenous variables key
open(8,file="endog_key.txt")
do ct=1,numendog
write(8,*) endog_key(ct,:)
end do !ct
close(8)

!write exogenous variables key
open(8,file="exog_key.txt")
do ct=1,numexog
write(8,*) exog_key(ct,:)
end do !ct
close(8)

!write endogenous variables pos key
open(8,file="endog_pos.txt")
do ct=1,numendog
write(8,*) endog_pos(ct,:)
end do !ct
close(8)

!write exogenous variables pos key
open(8,file="exog_pos.txt")
do ct=1,numexog
write(8,*) exog_pos(ct,:)
end do !ct
close(8)

!write agg prod series
open(8,file="asimpos.txt")
do t=1,numper
write(8,*) asimpos(t)
end do !t
close(8)

!write unc series
open(8,file="ssimpos.txt")
do t=1,numper
write(8,*) ssimpos(t)
end do !t
close(8)

!!!!END INITIALIZATIONS/GRID SETUPS

GEct=1
!!!!!VF ITERATION BLOCK

write(*,*) " "

!now, create all of the matrices useful later for setting up RHS returns

!output Ymat(z,a,k,l)
!$omp parallel private(zct,act,kct,lct,zval,aval,kval,lval)
!$omp do collapse(4)
do zct=1,znum
do act=1,anum
do kct=1,knum
do lct=1,lnum
    zval = z0(zct); aval=a0(act); kval = k0(kct); lval=l0(lct)
    Ymat(zct,act,kct,lct) = y(zval,aval,kval,lval,alpha,nu)
end do; !lct
end do; !kct
end do; !act
end do; !zct
!$omp end do nowait
!$omp end parallel

!investment Imat(k,k')
do kct=1,knum
do kprimect=1,knum
    kval = k0(kct); kprimeval = k0(kprimect)
    Imat(kct,kprimect) = kprimeval - (1.0-deltak) * kval
end do !kprimect
end do !kct

!capital AC ACkmat(z,a,k,l,k')
!$omp parallel private(zct,act,kct,lct,kprimect,zval,aval,kval,lval,kprimeval)
!$omp do collapse(5)
do zct=1,znum
do act=1,anum
do kct=1,knum
do lct=1,lnum
do kprimect=1,knum
    zval = z0(zct); aval=a0(act); kval = k0(kct); lval = l0(lct); kprimeval = k0(kprimect)
    ACkmat(zct,act,kct,lct,kprimect) = &
        ACk(zval,aval,kval,lval,alpha,nu,kprimeval,capirrev,capfix,deltak)
end do !kprimect
end do !lct
end do !kct
end do !act
end do !zct
!$omp end do nowait
!$omp end parallel

!labor AC AClmat(exogct,k,l,l_{-1},K)
!$omp parallel private(zct,act,kct,lct,lmin1ct,sct,fcstct,zval,aval,kval,lval,&
!$omp& lmin1val,wfcstval,smin1ct,exogct)
!$omp do collapse(5)
do exogct=1,numexog
do kct=1,knum
do lct=1,lnum
do lmin1ct=1,lnum
do fcstct=1,numfcst
    zct = exog_pos(exogct,1);
    act = exog_pos(exogct,2);
    sct = exog_pos(exogct,3);
    smin1ct = exog_pos(exogct,4);

    zval = z0(zct); aval=a0(act); kval = k0(kct); lval = l0(lct); lmin1val = l0(lmin1ct);
    wfcstval = wfcstmat(act,sct,smin1ct,fcstct);
    AClmat(exogct,kct,lct,lmin1ct,fcstct) = &
        ACl(zval,aval,kval,lval,alpha,nu,lmin1val,hirelin,firelin,labfix,deltan,wfcstval)

end do !fcstct
end do !lmin1ct
end do !lct
end do !kct
end do !exogct
!$omp end do nowait
!$omp end parallel

!wage bill WLmat(a,s,s_{-1},K,l)
!$omp parallel private(act,sct,smin1ct,fcstct,lct)
!$omp do collapse(5)
do act=1,anum
do sct=1,snum
do smin1ct=1,snum
do fcstct=1,numfcst
do lct=1,lnum
    WLmat(act,sct,smin1ct,fcstct,lct) = wfcstmat(act,sct,smin1ct,fcstct) * l0(lct)
end do !lct
end do !fcstct
end do !smin1ct
end do !sct
end do !act
!$omp end do
!$omp end parallel

write(*,"(A,F15.1,A)") "Finished setup of return matrices at ",omp_get_wtime()-start," seconds"


!if you want to do VFI, just read the policies/VF in from files
if (doVFI==1) then

write(*,*) "Doing VFI."


!initialize the VF and policies
if (GEct==1) then
    !in this case, initialize using a stupid guess
    Vold(:,:,:) = 0.0;
    V(:,:,:) = 0.0
    polmatold(:,:,:) = numendog/2
    polmat(:,:,:) = 0
    EVmat(:,:,:) = 0.0
else if (GEct>1) then
    !in this case, initialize using last GE iteration
    Vold = V
    polmatold = polmat
    polmat(:,:,:) = 0
    EVmat(:,:,:) = 0.0
end if

!this loop is the "VFI loop" itself, although technically this is policy iteration
do vfct=1,vfmaxit

    !!!Howard acceleration step
    do accelct=1,accelmaxit

        !set up the parallelization for the Howard acceleration construction
        !$omp parallel private(ct,endogct,exogct,fcstct,zct,act,sct,kct,lmin1ct,polstar,&
        !$omp& kprimect,lct,ind,weight,Vnextval,smin1ct)
        !$omp do
        do ct=1,numstates

            !extract states from loop index matrix
            endogct=loopind(ct,1); exogct=loopind(ct,2); fcstct=loopind(ct,3);

            !extract positions associated with these states
            zct = exog_pos(exogct,1); act = exog_pos(exogct,2); sct = exog_pos(exogct,3);
            smin1ct=exog_pos(exogct,4)
            kct = endog_pos(endogct,1); lmin1ct= endog_pos(endogct,2);

            !extract policy from polmatold
            polstar = polmatold(endogct,exogct,fcstct)

            !extract positions associated with the policies
            kprimect = endog_pos(polstar,1); lct = endog_pos(polstar,2)

            !initialize the RHS by forming the current period return
            V(endogct,exogct,fcstct) = pfcstmat(act,sct,smin1ct,fcstct) * ( Ymat(zct,act,kct,lct) &
                - ACkmat(zct,act,kct,lct,kprimect) - AClmat(exogct,kct,lct,lmin1ct,fcstct) &
                - Imat(kct,kprimect) - WLmat(act,sct,smin1ct,fcstct,lct) )
            !write(*,*) "V(endogct,exogct,fcstct) = ",V(endogct,exogct,fcstct)
            !find interpolation interval and weights for fcst agg capital
            ind = kbarfcstinds(act,sct,smin1ct,fcstct)
            weight = kbarfcstweights(act,sct,smin1ct,fcstct)
            !write(*,*) "ind = ",ind
            !write(*,*) "weight = ",weight
            !add discounted expected continuation value
            do exogprimect=1,numexog

                !what is the linearly interpolated continuation value?
                Vnextval = weight * Vold(polstar,exogprimect,ind+1) + &
                    (1.0 - weight)*Vold(polstar,exogprimect,ind)

                !now, add the discounted expected value associated with this realization
                V(endogct,exogct,fcstct) = V(endogct,exogct,fcstct) + &
                    beta * pr_mat(exogct,exogprimect) * Vnextval

            end do !exogprimect


        end do !ct
        !$omp end do nowait
        !$omp end parallel

        !initialize for the next round of the Howard VFI
        Vold = V

    end do !accelct


    !!!create EVmat for RHS optimization

    !set up the parallelization for the creation of the continuation value
    !$omp parallel private(ct,polct,exogct,fcstct,act,sct,ind,weight,exogprimect,&
    !$omp& Vnextval,smin1ct)
    !$omp do
    do ct=1,numstates

        !extract exogenous and policy states
        polct=loopind(ct,1); exogct=loopind(ct,2); fcstct = loopind(ct,3)

        !what are act and sct?
        act = exog_pos(exogct,2); sct = exog_pos(exogct,3); smin1ct = exog_pos(exogct,4)

        !find interpolation interval and weights for fcst agg capital
        ind = kbarfcstinds(act,sct,smin1ct,fcstct)
        weight = kbarfcstweights(act,sct,smin1ct,fcstct)

        !create expected value for the next period based on policy, exogenous transitions, and forecast kbar
        EVmat(exogct,polct,fcstct) = 0.0
        do exogprimect=1,numexog; !loop over exogenous realizations next period

            !what is the linearly interpolated continuation value?
            Vnextval = weight * Vold(polct,exogprimect,ind+1) + &
                    (1.0 - weight)*Vold(polct,exogprimect,ind)

            !add the weighted continuation value to the next period
            EVmat(exogct,polct,fcstct) = EVmat(exogct,polct,fcstct) + pr_mat(exogct,exogprimect) * Vnextval
        end do !exogprimect

    end do !ct
    !$omp end do nowait
    !$omp end parallel

    !!!RHS optimization step

    !set up the parallelization for the optimization of the Bellman RHS
    !$omp parallel private(ct,endogct,exogct,fcstct,zct,act,sct,kct,lmin1ct,&
    !$omp& polct,kprimect,lct,RHSvec,polstar,smin1ct)
    !$omp do
    do ct=1,numstates

       !extract states from loop index matrix
        endogct=loopind(ct,1); exogct=loopind(ct,2); fcstct=loopind(ct,3);

        !extract positions associated with these states
        zct = exog_pos(exogct,1); act = exog_pos(exogct,2); sct = exog_pos(exogct,3);
        smin1ct=exog_pos(exogct,4)
        kct = endog_pos(endogct,1); lmin1ct= endog_pos(endogct,2);

        !create RHS vector by looping over policies
        do polct=1,numendog

            !extract positions associated with the policies
            kprimect = endog_pos(polct,1); lct = endog_pos(polct,2)

            !initialize the RHS by forming the current period return
            RHSvec(polct) = pfcstmat(act,sct,smin1ct,fcstct) * ( Ymat(zct,act,kct,lct) &
                - ACkmat(zct,act,kct,lct,kprimect) - AClmat(exogct,kct,lct,lmin1ct,fcstct) &
                - Imat(kct,kprimect) - WLmat(act,sct,smin1ct,fcstct,lct) )

            !actually form the RHS by adding continuation value
            RHSvec(polct) = RHSvec(polct) + beta * EVmat(exogct,polct,fcstct)

        end do !polct

        !extract policies via selection of max payoff
        polstar = maxloc(RHSvec,1)
        polmat(endogct,exogct,fcstct) = polstar
        V(endogct,exogct,fcstct) = RHSvec(polstar)

    end do !ct
    !$omp end do nowait
    !$omp end parallel

    !compute errors and output the info
    vferror = maxval(abs(V-Vold))
    polerror = maxval(abs(polmat-polmatold))


    write(*,"(A,I3,A,F15.5,A)") "VF iter = ",vfct," in ",omp_get_wtime()-start," seconds"
    write(*,"(A,I3,A,F15.5)") "VF iter = ",vfct,", VF error = ",vferror
    write(*,"(A,I3,A,F15.5)") "VF iter = ",vfct,", policy error = ",polerror
    write(*,*) " "

    !exit criterion, then initialize for the next round of VFI
    !if ((vferror<vferrortol .and. polerror<vferrortol).or.(vferror<vferrortol)) exit
    if (polerror<vferrortol) exit
    Vold = V
    polmatold = polmat

end do !vfmaxit

!convert polmat to policies for capital and labor
do endogct=1,numendog
do exogct=1,numexog
do fcstct=1,numfcst
    !convert policies to actual capital and labor values
    kprime_key(endogct,exogct,fcstct) = endog_key(polmat(endogct,exogct,fcstct),1)
    kprime_pos(endogct,exogct,fcstct) = endog_pos(polmat(endogct,exogct,fcstct),1)

    !convert policies to capital and labor integer grid points
    lpol_key(endogct,exogct,fcstct) = endog_key(polmat(endogct,exogct,fcstct),2)
    lpol_pos(endogct,exogct,fcstct) = endog_pos(polmat(endogct,exogct,fcstct),2)
end do !fcstct
end do !exogct
end do !endogct

!display warning if hit the top of the grid or the bottom of the grid
!for either policy variable

write(*,*) " "


write(*,*) "Hit the top or bottom of grid in theory?"


if (maxval(kprime_pos)==knum) then
    write(*,*) "Hit top of capital grid."

end if

if (maxval(lpol_pos)==lnum) then
    write(*,*) "Hit top of labor grid."

end if

if (minval(kprime_pos)==1) then
    write(*,*) "Hit bottom of capital grid."

end if

if (minval(lpol_pos)==1) then
    write(*,*) "Hit bottom of labor grid."

end if

!if you don't want to do VFI, just read the policies/VF in from files
else if (doVFI==0) then

write(*,*) "Not doing VFI: reading value function and policies."


!read VF from file
open(8,file="V.txt")
do endogct=1,numendog
do exogct=1,numexog
do fcstct=1,numfcst
read(8,*) V(endogct,exogct,fcstct)
end do !fcstct
end do !exogct
end do !endogct
close(8)

!read policy matrix from file
open(8,file="polmat.txt")
do endogct=1,numendog
do exogct=1,numexog
do fcstct=1,numfcst
read(8,*) polmat(endogct,exogct,fcstct)
end do !fcstct
end do !exogct
end do !endogct
close(8)

!convert polmat to policies for capital and labor
do endogct=1,numendog
do exogct=1,numexog
do fcstct=1,numfcst
    !convert policies to actual capital and labor values
    kprime_key(endogct,exogct,fcstct) = endog_key(polmat(endogct,exogct,fcstct),1)
    kprime_pos(endogct,exogct,fcstct) = endog_pos(polmat(endogct,exogct,fcstct),1)

    !convert policies to capital and labor integer grid points
    lpol_key(endogct,exogct,fcstct) = endog_key(polmat(endogct,exogct,fcstct),2)
    lpol_pos(endogct,exogct,fcstct) = endog_pos(polmat(endogct,exogct,fcstct),2)
end do !fcstct
end do !exogct
end do !endogct

end if !doVFI

write(*,*) " "


!!!!!END VF ITERATION BLOCK

!!!!!SIMULATION BLOCK

write(*,*) "Doing unconditional simulation."

!first, initialize the distribution over endogenous variables
distzkl(:,:,:) = 0.0; !set everything to zero

if (distinit==1) then

    !if you read the initial distribution in
    open(8,file="distzklinit.txt")
    do zct=1,znum
    do endogct = 1,numendog
        read(8,*) distzkl(zct,endogct,1)
    end do !endogct
    end do !zct
    close(8)

else


    !chosen point to center the mass
    kct = 50
    lmin1ct=20
    endogct = (kct-1)*lnum + lmin1ct

    distzkl(:,(endogct-5):(endogct+5),1) = 1.0

end if !distinit

!now, round the dist
distzkl(:,:,1) = distzkl(:,:,1) / sum(distzkl(:,:,1))

!then, initialize the aggregate series to 0
Ysim(:) = 0.0; Ksim = Ysim; Lsim = Ysim; Isim = Ysim; Hsim = Ysim; ACksim = Ysim; AClsim = Ysim;
TOPksim = Ysim; BOTksim = Ysim; TOPlsim = Ysim; BOTlsim = Ysim; ADJksim = Ysim; ADJlsim = Ysim;
Csim = Ysim; psim =Ysim;

!start with aggregate capital guessed at some reasonable value
Ksim(1) = (kbarmin + kbarmax)/dble(2.0)

!initialize the endogenous state variable positions of firms
endogfirmpos(1,:) = numendog/2 !initialize the individual firms at the median endogenous grid point
!recall, the exogenous idio prod values are already determined

firmpos(:,:) = 0 !this array, (znum,numendog), is zero when no firms are at an index
do firmct=1,nfirms
	zct = zfirmpos(1,firmct)
	endogct = endogfirmpos(1,firmct)
	firmpos(zct,endogct) = firmpos(zct,endogct) + 1
end do !firmct


write(*,"(A,F15.1,A)") "Starting actual sim at ",omp_get_wtime()-start," seconds."


!now, loop over periods in the simulation, tracking the endogenous
!movement of weight around the grid
do t=1,numper-1

    !aggregate states
    if (t>1) then
        !in this case you don't get a segmentation fault when looking for s_{-1} index
        act = asimpos(t); sct=ssimpos(t); smin1ct = ssimpos(t-1)
    else if (t==1) then
        !initialize the very first one to low unc previously
        act = asimpos(t); sct=ssimpos(t); smin1ct = 1;
    end if
    aval = a0(act)

    !what is the next period forecast capital, in grid boundaries? index and weight?
    kbarfcstval = (kbar0(1)+kbar0(2))/dble(2.0)

    ind = 1
    weight=1.0

    Kfcstsim(t+1) = kbarfcstval

    !what is the current period forecast price?
    pfcstsim(t) = pval

	!set up the interpolation of the excess demand function e(p) = 1/p - C(p)
	ep0(:) = 0.0
	Cp0(:) = 0.0
	Yp0(:) = 0.0
	Ip0(:) = 0.0
	ACkp0(:) = 0.0
	AClp0(:) = 0.0
	Kbarprimep0(:) = 0.0
	Hp0(:) = 0.0
	Lp0(:) = 0.0
	polmatp_interp(:,:,:) = 0

		piter = 1
		pval = p0(piter) !loop over the vector of pre-stored p values
		wval = theta / pval; !what is the wage implied by this value of p?

        !this block computes the policies given the current pval & wval at each
        !point (z,kl_{-1})

        Yvalp = 0.0;
        Ivalp = 0.0;
        ACkvalp = 0.0;
        AClvalp = 0.0;
        Kbarprimevalp = 0.0;
        Hvalp = 0.0;
        Lvalp = 0.0

        !$omp parallel private(zct,endogct,exogct,kct,lmin1ct,RHSvec,polct,kprimect,&
        !$omp& lct,polstar) reduction(+:Yvalp,Ivalp,ACkvalp,AClvalp,Kbarprimevalp,Hvalp,&
        !$omp& Lvalp)
        !$omp do collapse(2)
        do zct=1,znum
        do endogct=1,numendog
            !restriction on which policies to compute - important for time
            if ((distzkl(zct,endogct,t)>disttol).or.(firmpos(zct,endogct)>0)) then
            !extract states from loop indexes
            exogct = (zct-1)*anum*snum*snum + (act-1)*snum*snum + (sct-1)*snum  + smin1ct

            !extract positions associated with these states
            kct = endog_pos(endogct,1);
            lmin1ct= endog_pos(endogct,2);

            !what is the new optimum value in index form and in capital and labor indexes?
            vmatp_interp(zct,endogct,piter) = V(endogct,exogct,1)/pval

            !what are the investment and labor for this position?
            lct = lpol_pos(endogct,exogct,1); !current labor, policy
            kprimect = kprime_pos(endogct,exogct,1); !next period capital, policy

            !what are the aggregates?
            Yvalp = Yvalp + distzkl(zct,endogct,t) * Ymat(zct,act,kct,lct)
            Ivalp = Ivalp + distzkl(zct,endogct,t) * Imat(kct,kprimect)
            ACkvalp = ACkvalp + distzkl(zct,endogct,t) * ACkmat(zct,act,kct,lct,kprimect)
            AClvalp = AClvalp + distzkl(zct,endogct,t) * ACl(z0(zct),a0(act),k0(kct),&
                l0(lct),alpha,nu,l0(lmin1ct),hirelin,firelin,labfix,deltan,wval)
            Kbarprimevalp = Kbarprimevalp + distzkl(zct,endogct,t) * k0(kprimect)
            Hvalp = Hvalp + distzkl(zct,endogct,t) * ( l0(lct) - (1.0-deltan) * l0(lmin1ct))
            Lvalp = Lvalp + distzkl(zct,endogct,t) * l0(lct)

            end if !disttol
        end do !endogct
        end do !zct
        !$omp end do nowait
        !$omp end parallel

        !what is implied consumption?
        Cvalp = Yvalp - Ivalp - ACkvalp - AClvalp

        !insert implied value of excess demand, and all other values into storage vectors
        ep0(piter) = 1.0/pval - Cvalp
        Cp0(piter) = Cvalp
       	Yp0(piter) = Yvalp
		Ip0(piter) = Ivalp
		ACkp0(piter) = ACkvalp
		AClp0(piter) = AClvalp
		Kbarprimep0(piter) = Kbarprimevalp
		Hp0(piter) = Hvalp
		Lp0(piter) = Lvalp


    !insert market-clearing price and other linearly interpolated stuff into sim series
    psim(t) = pval; !this is the most recently run p from the clearing algorithm
    Csim(t) = Cvalp; !this is already the linearly interpolated consumption C
    Ysim(t) = Yvalp; !linearly interpolated output Y
    Isim(t) = Ivalp; !linearly interpolated investment I
    ACksim(t) = ACkvalp; !linearly interpolated ACk
    AClsim(t) = ACLvalp; !linearly interpolated ACl
    Ksim(t+1) = Kbarprimevalp; !linearly interpolated K'
    Hsim(t) = Hvalp; !linearly interpolated hiring H
    Lsim(t) = Lvalp; !linearly interpolated labor input L

    !now that the market-clearing price is determined, move on to insert weight into the next period, according to the
    !linearly interpolated rule

    do zct=1,znum
    do endogct=1,numendog
        if (distzkl(zct,endogct,t)>disttol) then

         exogct = (zct-1)*anum*snum*snum + (act-1)*snum*snum + (sct-1)*snum  + smin1ct


        !based on the latest price, what is the policy here at pind?
        polstar = polmat(endogct,exogct,1)

        !insert distributional weight in appropriate slots next period
        distzkl(:,polstar,t+1) = distzkl(:,polstar,t+1) + &
        	pr_mat_z(zct,:,sct)*distzkl(zct,endogct,t)


        end if
    end do !endogct
    end do !zct

    !now, round to make sure that you're ending up with a distribution which makes sense each period
    distzkl(:,:,t+1) = distzkl(:,:,t+1)/sum(distzkl(:,:,t+1))


	!now, increment the firm simulation appropriately
	firmpos(:,:) = 0
	do firmct=1,nfirms

		!what is this firm's idio prod and current endogenous state?
		zct = zfirmpos(t,firmct)
		endogct = endogfirmpos(t,firmct)

		exogct = (zct-1)*anum*snum*snum + (act-1)*snum*snum + (sct-1)*snum  + smin1ct

		!what is the firms policy?
		polstar = polmat(endogct,exogct,1)
		endogfirmpos(t+1,firmct) = polstar

		!what is the firm's value?
		vfirm(t,firmct) = vmatp_interp(zct,endogct,1)

		!firm output
		lmin1ct = endog_pos(endogct,2)
		lct = endog_pos(polstar,2)
		kct = endog_pos(endogct,1)
		kprimect = endog_pos(polstar,1)
		yval = Ymat(zct,act,kct,lct)
		yfirm(t,firmct) = yval
		lfirm(t,firmct) = l0(lct)
		dfirm(t,firmct) = y(z0(zct),a0(act),k0(kct),l0(lct),alpha,nu)  &
		- ACk(z0(zct),a0(act),k0(kct),l0(lct),alpha,nu,k0(kprimect),capirrev,capfix,deltak) &
		- ACl(z0(zct),a0(act),k0(kct),l0(lct),alpha,nu,l0(lmin1ct),hirelin,firelin,labfix,deltan,wval) &
		- (k0(kprimect)-(1.0-deltak)*k0(kct)) &
		- wval * l0(lct)


		!now, increment the firm position array for next period
		zprimect = zfirmpos(t+1,firmct)
		firmpos(zprimect,polstar) = firmpos(zprimect,polstar) + 1

		!compute returns if possible
		if (t>1) then
			returnfirm(t,firmct) = log(vfirm(t,firmct)/(vfirm(t-1,firmct)-dfirm(t-1,firmct)))

		end if

		!if possible, compute within-firm return standard deviation
		if (t>4) then
			meanval =  (1.0/4.0)*(returnfirm(t,firmct)+returnfirm(t-1,firmct)+&
				returnfirm(t-2,firmct) + returnfirm(t-3,firmct))
			sdval = (1.0/4.0)*(returnfirm(t,firmct)**dble(2.0)+returnfirm(t-1,firmct)**dble(2.0)+&
				returnfirm(t-2,firmct)**dble(2.0)+ returnfirm(t-3,firmct)**dble(2.0))
			returnfirmsd(t,firmct) = sqrt(sdval - meanval**dble(2.0))
		end if

	end do !firmct


end do !t

!process firm level data to compute series of interest
FIRSTsim(:) = 0.0; SECONDsim = FIRSTsim; GDPsim = FIRSTsim; GROWTHsim = FIRSTsim;

!first, compute noise in the firm-level returns
returnmean = 0.0
returnstdev = 0.0
ct = 0
do t=2,numper-1
	do firmct=1,nfirmspub
		ct = ct+1
		returnmean = returnmean + returnfirm(t,firmct)
		returnstdev = returnstdev + returnfirm(t,firmct)**dble(2.0)
	end do !firmct
end do !t
returnmean = returnmean/dble(ct)
returnstdev = returnstdev/dble(ct)
returnstdev = sqrt(returnstdev-returnmean**dble(2.0))

!add noise to individual returns
returnfirmnoise = returnfirm + returnstdev * normfirmnoise


!now, recompute cross-sectional SD with noise
do t=4,numper-1
do firmct=1,nfirmspub
	meanval =  (1.0/4.0)*(returnfirmnoise(t,firmct)+returnfirmnoise(t-1,firmct)+&
				returnfirmnoise(t-2,firmct) + returnfirmnoise(t-3,firmct))
	sdval = (1.0/4.0)*(returnfirmnoise(t,firmct)**dble(2.0)+returnfirmnoise(t-1,firmct)**dble(2.0)+&
			returnfirmnoise(t-2,firmct)**dble(2.0)+ returnfirmnoise(t-3,firmct)**dble(2.0))
			returnfirmsdnoise(t,firmct) = sqrt(sdval - meanval**dble(2.0))
end do !firmct
end do !t

!compute the moments of interest
do t=2,numper-1

	do firmct=1,nfirms

		GDPsim(t) = GDPsim(t) + yfirm(t,firmct)

		if (t>1.and.firmct<=nfirmspub) then
			FIRSTsim(t) = FIRSTsim(t) + returnfirmnoise(t,firmct)
			SECONDsim(t) = SECONDsim(t) + returnfirmnoise(t,firmct)**dble(2.0)
		end if

	end do !firmct

	!normalize first and second moments of return distribution

	!THIS IS FRET
	FIRSTsim(t) = FIRSTsim(t)/dble(nfirmspub)

	!THIS IS SRETcs
	SECONDsim(t) = SECONDsim(t)/dble(nfirmspub)
	SECONDsim(t) = sqrt(SECONDsim(t)-FIRSTsim(t)**dble(2.0))

	if (t>1) then
		GROWTHsim(t) = dble(100.0)*log(GDPsim(t)/GDPsim(t-1))
	end if

end do !t

do t=5,numper-1

	!average of the CS returns for the year (micro concept)
	SECONDsimYR(t) = dble(0.25)*(SECONDsim(t)+SECONDsim(t-1)+SECONDsim(t-2)+SECONDsim(t-3))


end do !t

!now, do annual versions
GROWTHsimYR(:) = 0.0; FIRSTsimYR(:) = 0.0;

do t=5,numper-1

GROWTHsimYR(t) = dble(0.25)*(GROWTHsim(t)+GROWTHsim(t-1)+GROWTHsim(t-2)+GROWTHsim(t-3))

!THIS IS FRETann
FIRSTsimYR(t) = dble(0.25)*(FIRSTsim(t)+FIRSTsim(t-1)+FIRSTsim(t-2)+FIRSTsim(t-3))


	!square of the return for the year
	SECONDsim(t) = (FIRSTsim(t)-0.01)**dble(2.0)+&
		(FIRSTsim(t-1)-0.01)**dble(2.0)+&
		(FIRSTsim(t-2)-0.01)**dble(2.0)+&
		(FIRSTsim(t-3)-0.01)**dble(2.0)

end do !t

do t=5,numper-1

	!average of the average of returns for the year
	FIRSTsim(t) = FIRSTsimYR(t)

end do !t

write(*,*) " "
write(*,"(A,F15.1,A)") "Finished dist manipulation in ",omp_get_wtime()-start," seconds."


if (checkbounds==1) then

    write(*,*) " "
    write(*,*) "Checking for bounds of state space in simulation."



    TOPksim(:) = 0.0
    BOTksim(:) = 0.0
    TOPlsim(:) = 0.0
    BOTlsim(:) = 0.0

    do t = 1,numper-1
        do zct=1,znum
        do endogct=1,numendog
            kct = endog_pos(endogct,1); lmin1ct = endog_pos(endogct,2);
            if (kct==1) BOTksim(t) = BOTksim(t) + distzkl(zct,endogct,t)
            if (lmin1ct==1) BOTlsim(t) = BOTlsim(t) + distzkl(zct,endogct,t)
            if (kct==knum) TOPksim(t) = TOPksim(t) + distzkl(zct,endogct,t)
            if (lmin1ct==lnum) TOPlsim(t) = TOPlsim(t) + distzkl(zct,endogct,t)
        end do !endogct
        end do !zct

    end do !t

	write(*,*) "WEIGHT AT BOUNDARIES"
	write(*,*) maxval(BOTksim((numdiscard+1):(numper-1))),maxval(TOPksim((numdiscard+1):(numper-1))),&
		maxval(BOTlsim((numdiscard+1):(numper-1))),maxval(TOPlsim((numdiscard+1):(numper-1)))

	!compute the sum total of weight at the boundaries
	totweightval = maxval(BOTksim((numdiscard+1):(numper-1)))+maxval(TOPksim((numdiscard+1):(numper-1)))&
		+maxval(BOTlsim((numdiscard+1):(numper-1)))+maxval(TOPlsim((numdiscard+1):(numper-1)))

end if !checkbounds
write(*,*) " "
write(*,"(A,F15.1,A)") "Finished simulation in ",omp_get_wtime()-start," seconds."



!!!!!END SIMULATION BLOCK

!!!!WRITE THE PRELIM RESULTS BEFORE DETERMINING EXIT STATUS

write(*,*) "Writing output to .txt file for use in MATLAB."


write(*,*) " "



if (doVFI==1) then
!VF
open(8,file="V.txt")
do endogct=1,numendog
do exogct=1,numexog
do fcstct=1,numfcst
write(8,*) V(endogct,exogct,fcstct)
end do !fcstct
end do !exogct
end do !endogct
close(8)

!policy matrix
open(8,file="polmat.txt")
do endogct=1,numendog
do exogct=1,numexog
do fcstct=1,numfcst
write(8,*) polmat(endogct,exogct,fcstct)
end do !fcstct
end do !exogct
end do !endogct
close(8)

end if; !doVFI==1

!turn off the auxiliary file output for now
if (1==0) then
!output series
open(8,file="Ysim.txt")
do t=1,numper
write(8,*) Ysim(t)
end do !t
close(8)

!capital series
open(8,file="Ksim.txt")
do t=1,numper
write(8,*) Ksim(t)
end do !t
close(8)

!investment series
open(8,file="Isim.txt")
do t=1,numper
write(8,*) Isim(t)
end do !t
close(8)

!labor series
open(8,file="Lsim.txt")
do t=1,numper
write(8,*) Lsim(t)
end do !t
close(8)

!hiring series
open(8,file="Hsim.txt")
do t=1,numper
write(8,*) Hsim(t)
end do !t
close(8)

!consumption series
open(8,file="Csim.txt")
do t=1,numper
write(8,*) Csim(t)
end do !t
close(8)

!eqbm price series
open(8,file="psim.txt")
do t=1,numper
write(8,*) psim(t)
end do !t
close(8)


!capital AC series
open(8,file="ACksim.txt")
do t=1,numper
write(8,*) ACksim(t)
end do !t
close(8)

!labor AC series
open(8,file="AClsim.txt")
do t=1,numper
write(8,*) AClsim(t)
end do !t
close(8)

!top capital grid series
open(8,file="TOPksim.txt")
do t=1,numper
write(8,*) TOPksim(t)
end do !t
close(8)

!bottom capital grid series
open(8,file="BOTksim.txt")
do t=1,numper
write(8,*) BOTksim(t)
end do !t
close(8)

!top labor grid series
open(8,file="TOPlsim.txt")
do t=1,numper
write(8,*) TOPlsim(t)
end do !t
close(8)

!bottom labor grid series
open(8,file="BOTlsim.txt")
do t=1,numper
write(8,*) BOTlsim(t)
end do !t
close(8)


!forecast price series
open(8,file="pfcstsim.txt")
do t=1,numper
write(8,*) pfcstsim(t)
end do !t
close(8)

!forecast capital series
open(8,file="Kfcstsim.txt")
do t=1,numper
write(8,*) Kfcstsim(t)
end do !t
close(8)

!price iterations series
open(8,file="piterstoresim.txt")
do t=1,numper
write(8,*) piterstoresim(t)
end do !t
close(8)


!write firm value series
open(8,file="vfirm.txt")
do t=1,numper
do firmct=1,nfirms
write(8,*) vfirm(t,firmct)
end do !firmct
end do !t
close(8)

!write firm value added series
open(8,file="yfirm.txt")
do t=1,numper
do firmct=1,nfirms
write(8,*) yfirm(t,firmct)
end do !firmct
end do !t
close(8)

!write firm labor series
open(8,file="lfirm.txt")
do t=1,numper
do firmct=1,nfirms
write(8,*) lfirm(t,firmct)
end do !firmct
end do !t
close(8)



!first moment series
open(8,file="FIRSTsim.txt")
do t=1,numper
write(8,*) FIRSTsim(t)
end do !t
close(8)

!second moment series
open(8,file="SECONDsim.txt")
do t=1,numper
write(8,*) SECONDsim(t)
end do !t
close(8)

!GDP series
open(8,file="GDPsim.txt")
do t=1,numper
write(8,*) GDPsim(t)
end do !t
close(8)

!growth series
open(8,file="GROWTHsim.txt")
do t=1,numper
write(8,*) GROWTHsim(t)
end do !t
close(8)

!growth series year on year
open(8,file="GROWTHsimYR.txt")
do t=1,numper
write(8,*) GROWTHsimYR(t)
end do !t
close(8)

!return moments averaged over past year
open(8,file="FIRSTsimYR.txt")
do t=1,numper
write(8,*) FIRSTsimYR(t)
end do !t
close(8)

!return moments averaged over past year
open(8,file="SECONDsimYR.txt")
do t=1,numper
write(8,*) SECONDsimYR(t)
end do !t
close(8)


!write out the disaster files
open(8,file="NatDissim.txt")
open(9,file="PolShocksim.txt")
open(10,file="Revolutionsim.txt")
open(11,file="Terrorsim.txt")
do t=1,numper
write(8,*) DISASTERpos(t,1)
write(9,*) DISASTERpos(t,2)
write(10,*) DISASTERpos(t,3)
write(11,*) DISASTERpos(t,4)
end do !t
close(8)
close(9)
close(10)
close(11)

end if !end the switch for turning off auxiliary file output


outputfilestr = 'MATLABdata.csv'

!write out the file for processing via MATLAB
open(8,file=outputfilestr,recl=1000)
write(8,*) "x,growthsim,growthsimyr,fret,sretcs,fretann,srettsann,natdissim,polshocksim,revolutionsim,terrorsim,cflag,tsset_period"
do countryct=1,Ncountries
do t=1,Tper
ct = numdiscard + (countryct-1)*Tper + t
if (ct<numper) then
write(8,*) "1,",GROWTHsim(ct),&
	",",GROWTHsimYR(ct),&
	",",FIRSTsim(ct-1),&
	",",SECONDsim(ct-1),&
	",",FIRSTsimYR(ct-1),&
	",",SECONDsimYR(ct-1),&
	",",DISASTERpos(ct-1,1),&
	",",DISASTERpos(ct-1,2),&
	",",DISASTERpos(ct-1,3),&
	",",DISASTERpos(ct-1,4),&
	",",countryct,&
	",",t
end if
end do !t
end do !countryct
close(8)


write(*,*) "Calling MATLAB to run first and second-stage regressions."


write(*,*) " "


matlabstr = 'matlab -nodesktop -nosplash -r "FIRST_STAGE;'//'quit;"'

call system(matlabstr); !this actually calls MATLAB

!check to see if MATLAB has finished computing the moments, repeatedly
do ct=1,maxwaits

!give MATLAB waitsec seconds to complete the calculation
call sleep(waitsec)

write(*,"(A,I5,A)") "Checking ",ct," times for MATLAB completion."


!check to see if MATLAB's done
if (file_exists("FIRST_STAGE_OUT.txt")) exit

if (ct==maxwaits) then

write(*,"(A,I5,A)") "MATLAB failed to complete calculation in ",maxwaits," periods."


end if

end do !ct
write(*,*) "Done with MATLAB calculations."


write(*,*) " "


!now, insert the results into MATLAB_MOMS
MATLAB_MOMS(:) = 0.0
open(8,file="FIRST_STAGE_OUT.txt")
do ct=1,20
read(8,*) MATLAB_MOMS(ct)
end do !ct
close(8)


!data moment organization
!1 - nat dis to first, macro
!2 - pol shock to first, macro
!3 - revolution to first, macro
!4 - terror to first, macro
!5 - nat dis to second, macro
!6 - pol shock to second, macro
!7 - revolution to second, macro
!8 - terror to second, macro
!9 - second stage, first moment, macro
!10 - second stage, second moment, macro
!11 - nat dis to first, micro
!12 - pol shock to first, micro
!13 - revolution to first, micro
!14 - terror to first, micro
!15 - nat dis to second, micro
!16 - pol shock to second, micro
!17 - revolution to second, micro
!18 - terror to second, micro
!19 - second stage, first moment, micro
!20 - second stage, second moment, micro

!compute the scaled GMM objective
GMMobj = 0.0
write(*,*) "GMM Obj Relevant Moments (Data, Model, SE)"
do ct=1,20 !if both
!GMMobj = GMMobj + ((DATA_MOMS(ct) - MATLAB_MOMS(ct))/DATA_MOMS(ct))**dble(2.0)
!GMMobj = GMMobj + ((DATA_MOMS(ct) - MATLAB_MOMS(ct))/DATA_SE(ct))**dble(2.0)
if ((ct<=10).OR.(ct>14)) then
GMMobj = GMMobj + (DATA_MOMS(ct) - MATLAB_MOMS(ct))**dble(2.0)
end if
write(*,*) DATA_MOMS(ct),MATLAB_MOMS(ct)
end do !ct

write(*,*) " "

if (totweightval>=weighttol) GMMobj = 1000000000.0

!display the results
write(*,*) "Results of IV Regressions"

write(*,*) "MACRO FIRST STAGE, FIRST MOMENT"
write(*,*) "Data, Model"
do ct=1,4
write(*,*) DATA_MOMS(ct),MATLAB_MOMS(ct)
end do !ct

write(*,*) "MACRO FIRST STAGE, SECOND MOMENT"
write(*,*) "Data, Model"
do ct=5,8
write(*,*) DATA_MOMS(ct),MATLAB_MOMS(ct)
end do !ct

write(*,*) "MACRO SECOND STAGE"
write(*,*) "Data, Model"
do ct=9,10
write(*,*) DATA_MOMS(ct),MATLAB_MOMS(ct)
end do !ct

write(*,*) "MICRO FIRST STAGE, FIRST MOMENT"
write(*,*) "Data, Model"
do ct=11,14
write(*,*) DATA_MOMS(ct),MATLAB_MOMS(ct)
end do !ct

write(*,*) "MICRO FIRST STAGE, SECOND MOMENT"
write(*,*) "Data, Model"
do ct=15,18
write(*,*) DATA_MOMS(ct),MATLAB_MOMS(ct)
end do !ct

write(*,*) "MICRO SECOND STAGE"
write(*,*) "Data, Model"
do ct=19,20
write(*,*) DATA_MOMS(ct),MATLAB_MOMS(ct)
end do !ct

write(*,*) "GMM OBJ = ",GMMobj
fGMM = GMMobj

open(8,file="GMMmomstore.txt",position='append',recl=1000)
write(8,*) MATLAB_MOMS
close(8)

open(8,file="GMMobjstore.txt",position='append')
write(8,*) GMMobj
close(8)

open(8,file="GMMparamstore.txt",position='append',recl=1000)
write(8,*) DISASTERlev(1),DISASTERlev(2),DISASTERlev(3),DISASTERlev(4),&
	DISASTERuncprobs(1),DISASTERuncprobs(2),DISASTERuncprobs(3),DISASTERuncprobs(4)
close(8)

!call the clock to get the end time and closes log file
finish = omp_get_wtime()
write(*,"(A,F15.1,A)") "Finished completely in ",finish-start," seconds."

!do allocations for all variable-dimension arrays
deallocate(k0,l0,z0,a0,pr_mat_z,pr_mat_a,pr_mat_s,sigmaagrid,sigmazgrid,MATLABconstants,exog_key,exog_pos,endog_key,&
    endog_pos,V,pr_mat,EVmat,loopind,polmat,Vold,polmatold,RHSvec,kprime_key,lpol_key,kprime_pos,lpol_pos,asimpos,ssimpos,&
    distzkl,asimshocks,ssimshocks,Ysim,Ksim,Lsim,Isim,Hsim,ACksim,AClsim,TOPksim,BOTksim,TOPlsim,BOTlsim,ADJksim,&
		ADJlsim,psim,asimposIRF,ssimposIRF,asimshocksIRF,ssimshocksIRF,YsimIRF,KsimIRF,LsimIRF,IsimIRF,&
    HsimIRF,ACksimIRF,AClsimIRF,ADJksimIRF,ADJlsimIRF,kbar0,kbarfcstmat,pfcstmat,kfcstcoeff,pfcstcoeff,&
    wfcstmat,Ymat,Imat,ACkmat,AClmat,WLmat,kbarfcstweights,kbarfcstinds,EVmatp,polmatp,kprime_posp,&
    lpol_posp,Csim,Kfcstsim,pfcstsim,kcoeffstore,pcoeffstore,kfcstcoeffnew,pfcstcoeffnew,piterstoresim,&
    kfcsterrorstore,pfcsterrorstore,psimIRF,CsimIRF,KfcstsimIRF,pfcstsimIRF,kfcsterrorstoreIRF,&
		pfcsterrorstoreIRF,piterstoresimIRF,distzklbefore,distzklafter,pR2store,KR2store,pRMSEstore,KRMSEstore,&
    pSSE,KSSE,pSST,KSST,Kmean,pmean,KRMSEchangestore,pRMSEchangestore,KR2changestore,pR2changestore,&
		pDenHaanfcst,KDenHaanfcst,pmaxDenHaanstat,KmaxDenHaanstat,pavgDenHaanstat,KavgDenHaanstat,&
		pmaxDenHaanchangestore,KmaxDenHaanchangestore,pavgDenHaanchangestore,KavgDenHaanchangestore,&
    seedarray,p0,ep0,Cp0,Yp0,Ip0,ACkp0,AClp0,Kbarprimep0,Hp0,Lp0,polmatp_interp,vmatp_interp,&
    endogfirmpos,zfirmpos,yfirm,vfirm,zfirmshocks,pfirmshocks,FIRSTsim,SECONDsim,GDPsim,&
		GROWTHsim,firmpos,normfirmnoise,returnfirm,returnfirmnoise,GROWTHsimYR,FIRSTsimYR,SECONDsimYR,&
		DISASTERpos,DISASTERshocks,achgshocks,lfirm,schgshocks,returnfirmsd,dfirm,returnfirmsdnoise,&
    ascorrshocks,ushocks1,ushocks2,aggretnoise,normaggretnoise,ushocksagg1,ushocksagg2)

end function fGMM

subroutine unctauchen(transarray,grid,znum,zmin,zmax,snum,sigmagrid,rho,meanshift)
implicit none

!this routine discretizes the AR(1) productivity processes, subject to
!unc shocks

!note that this was double-checked - it duplicates fn_tranvarm2.m results

!input/output declarations
integer :: znum,snum
double precision :: zmin,zmax,rho
double precision :: transarray(znum,znum,snum),grid(znum),sigmagrid(snum)
double precision :: meanshift

!other declarations
integer :: zct,sct,zprimect
double precision :: log_grid(znum),standdev,gridinc

!create grid
call linspace(log_grid,log(zmin),log(zmax),znum)
grid = exp(log_grid)
gridinc = log_grid(2)-log_grid(1)

!loop over unc states
do sct=1,snum
    standdev = sigmagrid(sct)
    do zct=1,znum
        !middle intervals
        do zprimect=2,(znum-1)
            transarray(zct,zprimect,sct) = &
                normcdf(log_grid(zprimect)+gridinc/2.0,rho*log_grid(zct)+meanshift,standdev) - &
                normcdf(log_grid(zprimect)-gridinc/2.0,rho*log_grid(zct)+meanshift,standdev)
        end do !zprimect
        !first interval & last interval take remainder of mass
        transarray(zct,1,sct) = normcdf(log_grid(1)+gridinc/2.0,rho*log_grid(zct)+meanshift,standdev)
        transarray(zct,znum,sct) = 1.0 - normcdf(log_grid(znum)-gridinc/2.0,rho*log_grid(zct)+meanshift,standdev)
    end do !zct

    !impose that everything sums to 1 with rounding
    do zct=1,znum
        transarray(zct,:,sct) = transarray(zct,:,sct) / sum(transarray(zct,:,sct))
    end do !zct

end do !sct

end subroutine unctauchen

subroutine firmexogsim(zfirmshocks,zfirmpos,ssimpos,numper,nfirms,pr_mat_z,znum,zinit,snum)
implicit none

!this subroutine simulates the exogenous firm processes

!input/output declarations
integer :: zinit,znum,numper,nfirms,snum
double precision :: zfirmshocks(numper,nfirms),pr_mat_z(znum,znum,snum)
integer :: zfirmpos(numper,nfirms),ssimpos(numper)

!other declarations
integer :: t,firmct,ct,zct,zprimect,sct
double precision :: zsimgrid(znum)

!start with initialized zinit for the first period for all firms
zfirmpos(1,:) = zinit

!then, do remaining periods, iterating each firm separately
do t=1,numper-1

	sct = ssimpos(t) !what is today's uncertainty index?

	!which firm are you manipulating?
	do firmct=1,nfirms

		!what is their current idio prod index?
		zct = zfirmpos(t,firmct)

		!create vector of cumulative prob. thresholds
		zsimgrid(:) = 0.0
		zsimgrid(1) = pr_mat_z(zct,1,sct)
		do ct=2,znum
			zsimgrid(ct) = zsimgrid(ct-1) + pr_mat_z(zct,ct,sct)
		end do !zct
		zsimgrid(znum)=dble(1.0)

		!compare the firm shocks to prob thresholds, insert appropriate value
		if (zfirmshocks(t+1,firmct)<zsimgrid(1)) then
			zprimect=1
		else
			do ct=2,znum
				if (zsimgrid(ct-1)<=zfirmshocks(t+1,firmct).and.zsimgrid(ct)>zfirmshocks(t+1,firmct)) then
					zprimect=ct
				end if
			end do !ct
		end if

		!insert appropriate index for next period's productivity into the idio prod mat
		zfirmpos(t+1,firmct) = zprimect

	end do !firmct

end do !t



end subroutine firmexogsim

subroutine uncexogsimIRF(ssimshocksIRF,asimshocksIRF,ssimposIRF,asimposIRF,pr_mat_a,pr_mat_s,ainit,&
        lengthIRF,numsimIRF,shockperIRF,shocklengthIRF,anum,snum,sinit,singleshock)
implicit none

!this subroutine simulates the unc and agg prod processes in the IRF experiment

!input/output declarations
integer :: ainit,lengthIRF,numsimIRF,shockperIRF,shocklengthIRF,asimposIRF(lengthIRF,numsimIRF),&
    ssimposIRF(lengthIRF,numsimIRF),anum,snum,sinit,singleshock
double precision :: ssimshocksIRF(lengthIRF,numsimIRF),asimshocksIRF(lengthIRF,numsimIRF),&
    pr_mat_a(anum,anum,snum),pr_mat_s(snum,snum)

!other declarations
integer :: t,simct,perct,sprimect,ct,aprimect,dropflag
double precision :: ssimgrid(snum),asimgrid(anum)

!initialize the simulation for agg prod - unc is initialized below
asimposIRF(1,:) = ainit
ssimposIRF(1,:) = sinit

!!!START HERE IN THE LOW UNC FOLLOWED BY HIGH UNC CASE
if (singleshock==0) then
do simct=1,numsimIRF; !count IRF experiments
do t=2,lengthIRF; !count periods

    !do the initial, normal, simulation
    if (t>=shockperIRF.and.t<=(shockperIRF+shocklengthIRF)) then !this is the period of low uncertainty

        ssimposIRF(t,simct) = 1

    else if (t==(shockperIRF+shocklengthIRF+1)) then !this is the period of high uncertainty shock

        ssimposIRF(t,simct) = 2

    else !this means either before or after the shock: unc should evolve normally

        sprimect=0; !this will store new value

        !create vector of cumulative probability thresholds from transition mat
        ssimgrid(1) = pr_mat_s(ssimposIRF(t-1,simct),1)
        do ct=2,snum
            ssimgrid(ct) = ssimgrid(ct-1) + pr_mat_s(ssimposIRF(t-1,simct),ct)
        end do !ct

        !compare the uniform shock to the thresholds, and decide which bin
        if (ssimshocksIRF(t,simct)<ssimgrid(1)) then
            sprimect = 1
        else
            do ct=2,snum
                if (ssimgrid(ct-1)<=ssimshocksIRF(t,simct).and.ssimshocksIRF(t,simct)<ssimgrid(ct)) then
                    sprimect = ct
                end if
            end do !ct
        end if

        !fill in the simulated unc position
        ssimposIRF(t,simct) = sprimect

    end if
end do !t
end do !simct
end if !singleshock == 0



!!!START HERE IN THE SINGLE SHOCK CASE
if (singleshock==1) then
do simct=1,numsimIRF; !count IRF experiments
do t=2,lengthIRF; !count periods

    if (t==shockperIRF) then !this is the period of high uncertainty shock

        ssimposIRF(t,simct) = 2

    else !this means either before or after the shock: unc should evolve normally

        sprimect=0; !this will store new value

        !create vector of cumulative probability thresholds from transition mat
        ssimgrid(1) = pr_mat_s(ssimposIRF(t-1,simct),1)
        do ct=2,snum
            ssimgrid(ct) = ssimgrid(ct-1) + pr_mat_s(ssimposIRF(t-1,simct),ct)
        end do !ct

        !compare the uniform shock to the thresholds, and decide which bin
        if (ssimshocksIRF(t,simct)<ssimgrid(1)) then
            sprimect = 1
        else
            do ct=2,snum
                if (ssimgrid(ct-1)<=ssimshocksIRF(t,simct).and.ssimshocksIRF(t,simct)<ssimgrid(ct)) then
                    sprimect = ct
                end if
            end do !ct
        end if

        !fill in the simulated unc position
        ssimposIRF(t,simct) = sprimect

    end if
end do !t
end do !simct
end if !singleshock == 1

!now, simulate the agg prod process, according to the simulated unc process from above
do simct=1,numsimIRF
do t=2,lengthIRF; !note that this counter is over all the periods, not differentiating across simulations
    aprimect=0; !this will store the new value of agg prod

    !create vector of cumulative probability thresholds from transition mat
    !note that the third entry takes uncertainty into account
    asimgrid(1) = pr_mat_a(asimposIRF(t-1,simct),1,ssimposIRF(t-1,simct))
    do ct=2,anum
        asimgrid(ct) = asimgrid(ct-1) + pr_mat_a(asimposIRF(t-1,simct),ct,ssimposIRF(t-1,simct))
    end do !ct

    !compare the uniform shock to the thresholds, and decide which bin
    if (asimshocksIRF(t,simct)<asimgrid(1)) then
        aprimect = 1
    else
        do ct=2,anum
            if (asimgrid(ct-1)<=asimshocksIRF(t,simct).and.asimshocksIRF(t,simct)<asimgrid(ct)) then
                aprimect = ct
            end if
        end do !ct
    end if

    !fill in the simulated unc position
    asimposIRF(t,simct) = aprimect
end do !t
end do !simct

end subroutine uncexogsimIRF

subroutine uncexogsim(ssimshocks,asimshocks,DISASTERshocks,DISASTERprobs,ssimpos,asimpos,&
	pr_mat_s,pr_mat_a,ainit,sinit,snum,anum,numper,a0,achgshocks,DISASTERpos,sigmaa,schgshocks,&
	DISASTERuncprobs,DISASTERlev,ascorrshocks,firstsecondprob)
implicit none

!this subroutine simulates the uncertainty and aggregate productivity processes

!input/output declarations
integer :: numper,ainit,sinit,snum,anum,ssimpos(numper),asimpos(numper),&
	DISASTERpos(numper,4)
double precision :: ssimshocks(numper),asimshocks(numper),pr_mat_s(snum,snum),&
    pr_mat_a(anum,anum,snum),DISASTERshocks(numper,4),DISASTERprobs(4),a0(anum),&
    achgshocks(numper),sigmaa,schgshocks(numper),DISASTERuncprobs(4),&
    DISASTERlev(4),ascorrshocks(numper),firstsecondprob

!other declarations
integer :: t,ct,sprimect,aprimect,sct,act
double precision :: ssimgrid(snum),asimgrid(anum),ergtol,ergdistold(snum*anum),ergdist(snum*anum),ergerr,&
	ergdista(anum),ergmeana,totfirstmom,firstmomprob,secondmomprob

!first, compute joint ergodic distribution of (s,a)
ergtol=1.0e-6
ergdist(:) = 0.0
ergdistold(:) = 0.0
ergdistold(1) = 1.0
do t=1,1000

	do ct=1,snum*anum

		if (ct<=anum) then
			sct = 1
			act = ct
		else
			sct = 2
			act = ct - anum
		end if

		!tomorrow low unc?
		sprimect=1
		do aprimect=1,anum
			ergdist(aprimect) = ergdist(aprimect) + pr_mat_s(sct,sprimect)*pr_mat_a(act,aprimect,sct)*ergdistold(ct)
		end do !aprimect

		!tomorrow high unc?
		sprimect=2
		do aprimect=1,anum
			ergdist(aprimect+anum) = ergdist(aprimect+anum) + pr_mat_s(sct,sprimect)*pr_mat_a(act,aprimect,sct)*ergdistold(ct)
		end do !aprimect

	end do !ct


	ergerr = maxval(abs(ergdist-ergdistold))

	if (ergerr<ergtol) exit

	ergdistold = ergdist
	ergdist(:) = 0.0

end do !t

!now, extract ergodic distribution of a
ergdista(:) = 0.0
ergmeana = 0.0
do act=1,anum
	ergdista(act) = ergdist(act)+ergdist(act+anum)
	ergmeana = ergmeana + log(a0(act))*ergdista(act)
end do !act



!first, go through and determine whether disasters occur or not in each period
DISASTERpos(:,:) = 0
do t=1,numper
	do ct=1,4
		if (DISASTERshocks(t,ct)<=DISASTERprobs(ct)) DISASTERpos(t,ct)=1
	end do !ct
end do !t

!initialize the simulation for uncertainty and for agg prod
ssimpos(1) = sinit
asimpos(1) = ainit

!now simulate the uncertainty process
do t=2,numper
    sprimect=0; !this will store new value

    !create vector of cumulative probability thresholds from transition mat
    ssimgrid(1) = pr_mat_s(ssimpos(t-1),1)
    do ct=2,snum
        ssimgrid(ct) = ssimgrid(ct-1) + pr_mat_s(ssimpos(t-1),ct)
    end do !ct

    !compare the uniform shock to the thresholds, and decide which bin
    if (ssimshocks(t)<ssimgrid(1)) then
        sprimect = 1
    else
        do ct=2,snum
            if (ssimgrid(ct-1)<=ssimshocks(t).and.ssimshocks(t)<ssimgrid(ct)) then
                sprimect = ct
            end if
        end do !ct
    end if

    !fill in the simulated unc position
    ssimpos(t) = sprimect

    !now, if there is a disaster of type 2-4, i.e. if there is a coup, revolution, or terrorist
    !attack, uncertainty increases
    secondmomprob = 0.0
    if (DISASTERpos(t,1)==1) secondmomprob = secondmomprob + DISASTERuncprobs(1)
    if (DISASTERpos(t,2)==1) secondmomprob = secondmomprob + DISASTERuncprobs(2)
    if (DISASTERpos(t,3)==1) secondmomprob = secondmomprob + DISASTERuncprobs(3)
    if (DISASTERpos(t,4)==1) secondmomprob = secondmomprob + DISASTERuncprobs(4)

	if (schgshocks(t)<=secondmomprob) ssimpos(t) = 2

end do !t

!given the uncertainty process, simulate the aggregate productivity process
do t=2,numper
    aprimect=0; !this will store the new value of agg prod

    !create vector of cumulative probability thresholds from transition mat
    !note that the third entry takes uncertainty into account
    asimgrid(1) = pr_mat_a(asimpos(t-1),1,ssimpos(t-1))
    do ct=2,anum
        asimgrid(ct) = asimgrid(ct-1) + pr_mat_a(asimpos(t-1),ct,ssimpos(t-1))
    end do !ct

    !compare the uniform shock to the thresholds, and decide which bin
    if (asimshocks(t)<asimgrid(1)) then
        aprimect = 1
    else
        do ct=2,anum
            if (asimgrid(ct-1)<=asimshocks(t).and.asimshocks(t)<asimgrid(ct)) then
                aprimect = ct
            end if
        end do !ct
    end if

    !fill in the simulated unc position
    asimpos(t) = aprimect

	!now, is there a disaster of any sort?
	if (sum(DISASTERpos(t,:))>0) then

		!if so, figure out total first-moment impact
		totfirstmom = 0.0
		if (DISASTERpos(t,1)>0) totfirstmom = totfirstmom + sigmaa*DISASTERlev(1)
		if (DISASTERpos(t,2)>0) totfirstmom = totfirstmom + sigmaa*DISASTERlev(2)
		if (DISASTERpos(t,3)>0) totfirstmom = totfirstmom + sigmaa*DISASTERlev(3)
		if (DISASTERpos(t,4)>0) totfirstmom = totfirstmom + sigmaa*DISASTERlev(4)

		if (totfirstmom>=0.0) then

			firstmomprob = totfirstmom/(log(a0(anum))-ergmeana)
			if (achgshocks(t)<=firstmomprob) asimpos(t) = anum

		else

			firstmomprob = totfirstmom/(log(a0(1))-ergmeana)
			if (achgshocks(t)<=firstmomprob) asimpos(t) = 1

		end if


	end if


	!now, add correlation of second and first moments
	if ((ascorrshocks(t)<=firstsecondprob).and.(ssimpos(t)==2)) asimpos(t) = 1



end do !t
end subroutine uncexogsim

double precision function y(zval,aval,kval,lval,alpha,nu)
implicit none

!this function evaluates output given parameters and states

double precision :: zval,aval,kval,lval,alpha,nu

y = zval * aval * (kval ** alpha) * (lval ** nu)

end function

double precision function ACk(zval,aval,kval,lval,alpha,nu,kprimeval,capirrev,capfix,deltak)
implicit none

!this function evaluates capital AC given parameters and states

!input parameters
double precision :: zval,aval,kval,lval,alpha,nu,kprimeval,capirrev,capfix,deltak

!other parameters
double precision :: yval,changetol,ival

changetol = 1.0e-10
ival = kprimeval - (1.0-deltak) * kval

!start the AC at 0
ACk = 0.0

!take out the partial irreversibility costs
if (ival<-changetol) then
    ACk = ACk - ival * capirrev
end if

!take out the fixed disruption costs
if (abs(ival)>changetol) then
    yval = y(zval,aval,kval,lval,alpha,nu)
    ACk = ACk + yval * capfix
end if

end function ACk

double precision function ACl(zval,aval,kval,lval,alpha,nu,lmin1val,hirelin,firelin,labfix,deltan,wval)
implicit none

!this function evaluates labor AC given parameters and states

!input parameters
double precision :: zval,aval,kval,lval,alpha,nu,lmin1val,hirelin,firelin,labfix,deltan,wval

!other parameters
double precision :: hval,changetol,yval

changetol = 1.0e-10

hval = lval - (1.0-deltan)*lmin1val

!take out the fixed and linear costs if you changed costs
ACl = 0.0
if (abs(hval)>changetol) then
    yval = y(zval,aval,kval,lval,alpha,nu)
    ACl = ACl + labfix * yval; !the fixed costs
    if (hval<-changetol) ACl = ACl - hval * wval * firelin; !the linear firing costs
    if (hval>changetol) ACl = ACl + hval * wval * hirelin; !the linear hiring costs
end if
end function ACl


subroutine boxmuller(normshocks,M,N,ushocks1,ushocks2)

!input/output declarations
integer :: M,N
double precision :: normshocks(M,N),ushocks1(M,N),ushocks2(M,N)

integer :: ct1,ct2
double precision :: pi

!define pi
pi = 4.0*atan(1.0)

do ct1=1,M
do ct2=1,N

	normshocks(ct1,ct2) = sqrt(dble(-2.0)*log(ushocks1(ct1,ct2))) * cos(dble(2.0)*pi*ushocks2(ct1,ct2))

end do !ct2
end do !ct1

end subroutine boxmuller

logical function file_exists(filename)
implicit none

character (len=*) :: filename

inquire(file=filename,exist=file_exists)

end function file_exists


end program VOL_GROWTH_wrapper
