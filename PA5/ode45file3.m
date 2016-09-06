% Numerical Mathematics and Computing, Fifth Edition
% Ward Cheney & David Kincaid
% Brooks/Cole Publ. Co.
% (c) 2003
%
% file:  ode45file1.m
%
% This M-file is invoked by the file 'RKF_ode45.m'.

function xprime = ode45file3(t,x) 
r = .1;
g = 32.1;
% A = 512*pi/x;
% A = pi*(8-7.9/8*(x-8))^2;
A = pi*x^2;
xprime = -.6*pi*r^2*sqrt(2*g)*sqrt(x)/A;
% xprime_prime