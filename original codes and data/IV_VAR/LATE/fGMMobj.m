%%%%%%%
% fGMMobj.m
%
% This MATLAB function evaluates the GMM objective for the
% IV-VAR estimation in
% "Using Disasters to Estimate the Impact of Uncertainty,"
% by Scott R. Baker, Nicholas Bloom, and Stephen J. Terry.
%
% Initially run on MATLAB R2021a
%
% Questions to Stephen Terry
% stephenjamesterry@gmail.com
%%%%%%%
function OUT = fGMMobj(x)
global NX ND Nparams Nmoms MOMvec extraoutput

B = zeros(NX,NX);
Dcoeff = zeros(ND,2);

%extract the input values of the B matrix
B(1,1) = x(1); 
B(2,1) = x(2);
B(3,1) = x(3);
B(1,2) = x(4); 
B(2,2) = x(5);
B(3,2) = x(6);
B(1,3) = x(7); 
B(2,3) = x(8);
B(3,3) = x(9);

%extract the input values of the IV first stage matrix
Dcoeff(1,1) = x(10);
Dcoeff(2,1) = x(11);
Dcoeff(3,1) = x(12);
Dcoeff(4,1) = x(13);
Dcoeff(1,2) = x(14);
Dcoeff(2,2) = x(15);
Dcoeff(3,2) = x(16);
Dcoeff(4,2) = x(17);

LAMBDAimplied = [1 0 0; ...
    0 sum(Dcoeff(:,1).^2)+1 sum(Dcoeff(:,1).*Dcoeff(:,2));...
    0 sum(Dcoeff(:,1).*Dcoeff(:,2)) sum(Dcoeff(:,2).^2)+1];

%form the implied covariance matrix of the reduced-form residuals
OMEGAimplied = B * LAMBDAimplied * B';


%form the implied covariance moments of the reduced-form residsuals and the
%instruments
EDetaimplied =zeros(NX,ND);
for IVct=1:ND;
    EDetaimplied(:,IVct) = B*[0;Dcoeff(IVct,1); Dcoeff(IVct,2)];
end;


%compute implied moments
MOMimplied = zeros(Nmoms,1);
MOMimplied(1) = OMEGAimplied(1,1);
MOMimplied(2) = OMEGAimplied(2,2);
MOMimplied(3) = OMEGAimplied(3,3);
MOMimplied(4) = OMEGAimplied(1,2);
MOMimplied(5) = OMEGAimplied(1,3);
MOMimplied(6) = OMEGAimplied(2,3);
MOMimplied(7) = EDetaimplied(1,1);
MOMimplied(8) = EDetaimplied(2,1);
MOMimplied(9) = EDetaimplied(3,1);
MOMimplied(10) = EDetaimplied(1,2);
MOMimplied(11) = EDetaimplied(2,2);
MOMimplied(12) = EDetaimplied(3,2);
MOMimplied(13) = EDetaimplied(1,3);
MOMimplied(14) = EDetaimplied(2,3);
MOMimplied(15) = EDetaimplied(3,3);
MOMimplied(16) = EDetaimplied(1,4);
MOMimplied(17) = EDetaimplied(2,4);
MOMimplied(18) = EDetaimplied(3,4);



LHS = MOMvec;
RHS = MOMimplied;
GMMerrvec = LHS-RHS;
if (extraoutput==1);
    OUT = MOMimplied;
else
OUT = sum(GMMerrvec.^2);
end;

