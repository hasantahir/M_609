% MATH 609 - Programming Assignment 1
% Computational Exercise 3
% This program solves a linear system with a Hilbert matrix
close all;clc;clear
%% Matrix Definition
n = 20;
A = zeros(n,n);
x = zeros(n,1);
b = ones(n,1);
for i = 1 : n
    for j = 1 : n
        A(i,j) = 1/(i+j+1); % Fill up the Hilbert Matrix
    end
end
spyc(A*1000,30)
%% LU Factorization
[L,U] = lu_decompose(A);
%% Now follow Algorithm from "Numerical Methhods and Computing"
% Kincaid, Cheney, 6th ed, chapter 8, pg. 301:
z = zeros(n,1);
temp_z = 0;
temp_u = 0;
%% Evaluating Lc = b
z(1) = b(1);
for i = 2 : n
    temp_z = b(i);
    for j = 1 : i-1
        temp_z = temp_z - L(i,j)*z(j);
    end
    
    z(i) = temp_z;
end
%% Evaluating Ux = c
x(n) = z(n)/U(n,n);
for i = n-1 : -1 : 1
    temp_u = z(i);
    for j = i+1 : n
        temp_u = temp_u - U(i,j)*x(j);
    end
    
    x(i) = temp_u/U(i,i);
end
plot(x)
set(gcf,'Color','white');
plotTickLatex2D
% Create title       
matFileName = sprintf('math609_pa1_comp_example_3_n_%d', n);
saveas(gcf,[matFileName,'.eps'],'epsc')

