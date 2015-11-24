% Main program to evaluate cubic spline
% Hasan Tahir Abbas
% 01/11/2015
% Math 609 : Programming Assignment 4
% 
%
% This program reads data files and interpolates using cubic spline
% interpolation
clear;close all
%% Read data from the input file
filename = 'Car_travel.txt';
file_read = importdata(filename);
t = file_read(:,1); % the first column has the t vector of length N+1
y = file_read(:,2); % the second column has the y vector of length N+1
N = length(t); % Number of knots
z = spline3_coeff(N,t,y,'free');
x = 10;
S = zeros(1,length(x));
for i = 1 : length(x)
    S(i) = spline3_eval(N,t,y,z,x(i));
end
close all
plot(t,y,'o')
hold on
plot(x,S,'x','MarkerSize',10)