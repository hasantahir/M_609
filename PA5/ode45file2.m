% Numerical Mathematics and Computing, Fifth Edition
% Ward Cheney & David Kincaid
% Brooks/Cole Publ. Co.
% (c) 2003
%
% file:  ode45file1.m
%
% This M-file is invoked by the file 'RKF_ode45.m'.

function xprime = ode45file2(t,x) 

C = .3;
R = 1.4;
L = 1.7;
E = exp(-.06*pi*t).*sin(2*t-pi);      % range
E1 = (3*pi*sin(2*t)*exp(-(3*pi*t)/50))/50 - 2*cos(2*t)*exp(-(3*pi*t)/50);   % first derivative
E2 = 4*sin(2*t)*exp(-(3*pi*t)/50) + (6*pi*cos(2*t)*exp(-(3*pi*t)/50))/25 ...
    - (9*pi^2*sin(2*t)*exp(-(3*pi*t)/50))/2500;  % second derivative
xprime = C*E2 + 1/R*E1 + 1/L*E;


% (3*pi*sin(2*t)*exp(-(3*pi*t)/50))/50 - 2*cos(2*t)*exp(-(3*pi*t)/50)
% 
% 4*sin(2*t)*exp(-(3*pi*t)/50) + (6*pi*cos(2*t)*exp(-(3*pi*t)/50))/25...
% - (9*pi^2*sin(2*t)*exp(-(3*pi*t)/50))/2500