function [s t] = skewTestAC(x)
%function [s t] = skewTestAC(x)
%This function computes the coefficient of skewness (s) of x and gives a standard
%normally distributed test statistic (t) to test that s = 0.
%It follows Bai and Ng (JBES, 2005).
%Author: Alisdair McKay
%Date: 12/4/06

T = length(x);

s = skewness(x);


mu = mean(x);
temp = [(x-mu).^3  x-mu];

%gamma = parzKern(temp)
[gamma b] = lr_var(temp);
clear b;


a = [1 -3*var(x)];

t = s * sqrt(T/((a*gamma*a')/(std(x)^6)));
