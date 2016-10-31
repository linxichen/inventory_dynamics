function S=parzKern(X,nlags)
%function S=bartKern(X,nlags)
%This function estimates the spectral density of the process in the rows of
%X at freq. zero using the Parzen kernel.  nlags is the number of lags to
%use (default to T^(2/5), which yields fastest convergence according to
%Andrews(1991).

[T,k] = size(X);

if nargin < 2
    nlags = round(T^(2/5));
end

S = zeros(k);

for lag = 0:nlags
  rho = X(1:T-lag,:)'*X(1+lag:T,:)/(T-k);  
  if lag >= 1, rho = rho + rho'; end
  qx = lag/nlags;
  if qx <= 0.5
      wt = 1-6*qx^2+6*(abs(qx))^3;
  else
      wt = 2*(1-abs(qx))^3;
  end
  S = S + wt*rho;    
end

