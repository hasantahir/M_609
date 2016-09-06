% Numerical Mathematics and Computing, Fifth Edition
% Ward Cheney & David Kincaid
% Brooks/Cole Publ. Co.
% (c) 2003
%
% file:  ode45file1.m
%
% This M-file is invoked by the file 'RKF_ode45.m'.

function xprime = ode45file1(t,x) 
xprime = 2./t*x + t.^2*exp(t);
% xprime_prime