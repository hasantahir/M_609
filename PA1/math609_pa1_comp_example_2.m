% MATH 609 - Programming Assignment 1
% Computational Exercise 2
% This program solves a linear system with a Hilbert matrix
close all;clc;clear
%%
n = 40; % Choose n from 20 or 40
nd = n-2;
N = nd^2;
h = 1/(n+1);
A = zeros(N,N);
b = h^2*ones(n^2,1);
x = zeros(n^2,1);

%% Matrix Definition for 5-point finite difference equation
%% From Numerical Recipes
    for i = 2 : N -1
        for j = 2 : N -1
            if i == j
                A(i,j) = 4+h^2;
            elseif j == i-1 || j == i+1
                A(i,j) = -1;
            elseif j == i-n || j == i+n
                A(i,j) = -1;
            end
        end
    end
    for i = n:n:n^2-n
        A(i,i+1) = 0;
        A(i+1,i) = 0;
    end
        
spyc(A,30)
%% LU Factorization
% A = (A(2:n^2+1,2:n^2+1));
size(A)
[L,U] = lu_decompose(A);
%% Now follow Algorithm from "Numerical Methhods and Computing
% Kincaid, Cheney, 6th ed, chapter 8, pg. 301:
z = zeros(n^2-1,1);
temp_z = 0;
temp_u = 0;
%% Evaluating Lc = b
z(1) = b(1);
for i = 2 : n^2-1
    temp_z = b(i);
    for j = 1 : i-1
        temp_z = temp_z - L(i,j)*z(j);
    end
    
    z(i) = temp_z;
end
%% Evaluating Ux = c
x(n^2-1) = z(n^2-1)/U(n^2-1,n^2-1);
for i = n^2-1 : -1 : 1
    temp_u = z(i);
    for j = i+1 : n^2-1
        temp_u = temp_u - U(i,j)*x(j);
    end
    
    x(i) = temp_u/U(i,i);
end

x = reshape(x,[n,n]);

plot(1:n,(x(5,:)),'ko',1:n,(x(5,:)),'k-','LineWidth',1.2); %Piece-wise linear interpolation of x
set(gcf,'Color','white');
plotTickLatex2D
% Create title       
matFileName = sprintf('math609_pa1_comp_example_2_n_%d', n);
saveas(gcf,[matFileName,'.eps'],'epsc')