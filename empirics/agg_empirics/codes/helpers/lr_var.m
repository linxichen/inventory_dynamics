
% Bruce E. Hansen
% Department of Economics
% Social Science Building
% University of Wisconsin
% Madison, WI 53706-1393
% bhansen@ssc.wisc.edu
% http://www.ssc.wisc.edu/~bhansen/


% This file contains four procedures:
% UR_REG, LR_VAR, UR_CRITS, UR_ADF.



% (1) UR_REG produces the coefficient and studentized statistics
% for unit roots outlined in B. Hansen "Rethinking the univariate approach
% to unit root testing: Using covariates to increase power" 
% Econometric Theory (1995).


% The format is:


% {tstat,crits,rho2,sig2} = UR_REG(y,x,p,q,k1,k2);


% The inputs are:


% y   -  dependent variable, Tx1.
% x   -  right-hand-side variables, Txk.
% p   -  order of deterministic time trend.
%        if p=-1, then none included
%        if p=0,  then constant included
%        if p =1, then constant and time trend included.
% q   -  number of lagged dy terms
% k1  -  number of lagged x terms
% k2  -  number of lead x terms


% The outputs are:


% tstat   -  studentized coefficient test
% crits -  asymptotic critical values (1%, 5%, 10%) for tstat
% rho2    -  estimated "Rho-Squared"
% sig2    -  estimated error variance


% Notes:


% If k1=k2=0, then the variable x appears on the right-hand-side of the
% regression.  This is useful if the user wishes to directly control the
% context of the included covariates.


% The variable "x" needs to be stationary.  In most cases, it will be the
% first-difference of some levels variable.  (In the notation of the paper,
% it should be "dx").



% (2) LR_VAR  calculates a "Long-Run Variance" for a vector-valued
% series, using the suggestions of Andrews (1991).  The format is:


% {omega,bandw} = LR_VAR(e);


% where e is a Txk matrix of observations, omega is a kxk covariance
% matrix, and bandw is the bandwidth used in estimation of omega.



% Both UR_REG and LR_VAR depend on the following set of global variables.


% kernel_   -  indicator for kernel method
%              if kernel_ = 1, then Parzen kernel is used
%              if kernel_ = 2, then Quadratic Spectral kernel is used
%              if kernel_ = 3, then Bartlett (Newey-West) kernel is used
% band_     -  bandwidth used to calculate long-run variance matrix
%              if band_ = 0, then the bandwidth is calculated using
%              the suggestion of Andrews (1991)
% white_    -  pre-whiten error indicator
%              if white_ = 1, then the errors are pre-whitened with a VAR(1)
%              before kernel is applied.
% urprint_  -  print to screen supressor
%              if urprint_ = 0, then results (for UR_REG and UR_ADF) are printed
%              to screen
%              if urprint_ = 1, then no printing to screen
%             (useful for Monte Carlo runs)


% Global Defaults are set below
% Since they are globals, they may be changed in a program if used.



% (3)  UR_CRITS gives asymptotic critical values for the unit root tests.
% The format is


% crits = ur_crits(rho2,p);


% The inputs are:


% rho2  - R-squared measure, as given in UR_REG output
% p     -  order of deterministic time trend, as given in UR_REG input
%          if p=-1, then none included
%          if p=0,  then constant included
%          if p =1, then constant and time trend included.


% The outputs are:


% crits - asymptotic 1%, 5%, 10% critical values for tstat



% (4)  UR_ADF


% Format:


% {tstat,sig2} = UR_ADF(y,p,q);


% The inputs are:


% y   -  dependent variable, Tx1.
% p   -  order of deterministic time trend.
%        if p=-1, then none included
%        if p=0,  then constant included
%        if p =1, then constant and time trend included.
% q   -  number of lagged dy terms


% The outputs are:


% tstat  -  studentized coefficient test
% sig2   -  estimated error variance


% Print-to-screen is determined by the global variable _urprint,
% documented above.


% 


% PROCEDURE CODE: %



function [omega,bandw]=lr_var(u)
global kernel_;
global band_;
global urprint_;
global white_;
kernel_  = 1;
band_    = 0;
urprint_ = 0;
white_   = 1;


[tu,p]=size(u);
if white_==1
    te=tu-1;
    au=(u(2:tu,:)'/u(1:te,:)')';
    e=u(2:tu,:)-u(1:te,:)*au;
else
    e=u;
    te=tu;
end;


if band_==0
    eb=e(1:te-1,:);
    ef=e(2:te,:);
    ae=sum(eb.*ef)'./sum(eb.^2)';
    ee=ef-eb.*(ones(length(eb(:,1)),1)*ae');
    se=mean(ee.^2)';
    temp=(1-ae).^2;
    ad=sum((se./((1-ae).^2)).^2)';
    aese=ae.*se;
    a1=4*sum((aese./((((1-ae).^3).*(1+ae))*ones(1,length(aese(1,:))))).^2)'/ad;
    a2=4*sum((aese./(((1-ae).^4)*ones(1,length(aese(1,:))))).^2)'/ad;
    (aese./(((1-ae).^4)*ones(1,length(aese(1,:))))).^2;
    
    if kernel_==2
        bandw = 1.3221*((a2*te)^.2); % Quadratic Spectral %
    elseif kernel_==1
        bandw = 2.6614*((a2*te)^.2); % Parzen %
        if bandw>(te-2)
            bandw=te-2;
        end;
    elseif kernel_==3;               % Bartlett %
        bandw=1.1447*((a1*te)^.333);
        if bandw>(te-2)
            bandw=te-2;
        end;
    end;
else
    bandw=band_;
end;


% Estimate Covariances %        
if kernel_==2  % Quadratic Spectral Kernel %
    tm=te-1;
    jb=((1:tm)/bandw')';
    jband = jb*1.2*pi;
    kern = ((sin(jband)./jband - cos(jband))./(jband.^2)).*3;
elseif kernel_==1 % Parzen kernel %
    tm=floor(bandw);
    if tm>0
        jb=((1:tm)/bandw')';
        kern=(1-(jb.^2)*6+(jb.^3)*6).*(jb <=.5);
        kern=kern + ((1-jb).^3).*(jb>.5)*2;
    end;
elseif kernel_==3 % Bartlett kernel %
    tm=floor(bandw);
    if tm>0
        kern=1-((1:tm)/bandw')';
    end;
end;


lam=zeros(p,p);
j=1;
while j<=tm
    kj=kern(j);
    lam=lam+(e(1:te-j,:)'*e(1+j:te,:))*kj;
    j=j+1;
end;
omega=(e'*e+lam+lam')/te;


if white_==1
    eau=inv(eye(p)-au);
    omega=eau'*omega*eau;
end;
