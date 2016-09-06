% Numerical Mathematics and Computing, Fifth Edition
% Ward Cheney & David Kincaid
% Brooks/Cole Publ. Co.
% (c) 2003
%
% file:  ode45file1.m
%
% This M-file is invoked by the file 'RKF_ode45.m'.

function xprime = ode45file4(t,x) 
k = 6.22e-19;
n1 = 1000;
n2 = 1000;
n3 = 1500;
xprime = k*(n1-.5*x)^2*(n2-.5*x)^2*(n3-.75*x)^3;
% xprime_prime