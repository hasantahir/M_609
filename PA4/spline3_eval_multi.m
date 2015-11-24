function [S,S1,S2] = spline3_eval_multi(N,t,y,z,X)
% function z = spline3_eval(filename,x);
% Based on the pseudocode spline3_coeff presented in Kincaid and Cheney,
% " Numerical Mathematics and Computing"
%
% This function evaluates the cubic slpine at x 
%  INPUT:
%    filename  - data file containing knot and interpolation points as two
%                column vectors
%    X - point at which S(x) is to be evaluated
%  OUTPUT:
%     S - value of spline at x

%% Read data from the input file
% file_read = importdata(filename);
% t = file_read(:,1) % the first column has the t vector of length N+1
% y = file_read(:,2) % the second column has the y vector of length N+1
%% determine lengths of the intervals
n = N - 1;

for  i = n : -1 : 2
    
    if ( X - t(i)) >= 0
        break
    end
end

h = t(i+1) - t(i);
tmp = (z(i)/2) + (X - t(i))*(z(i+1) - z(i))/(6*h);
tmp = -(h/6)*(z(i+1) + 2*z(i)) + (y(i+1) - y(i))/h + (X - t(i))*tmp;

%% spline function
S = y(i) + (X - t(i))*tmp;

%% first derivative
S1 = z(i+1)/(2*h)*(X - t(i))^2 - z(i)/(2*h)*(t(i+1) - X)^2 ...
    + y(i+1)/h - h*z(i+1)/6 - y(i)/h + h*z(i)/6;

%% Second derivative
S2 = z(i+1)/h*(X - t(i)) + z(i)/h*(t(i+1) - X);

end
    