% MATH 609 - Programming Assignment 1
% Computational Exercise 1
% This program solves a linear system with a Hilbert matrix
close all;clc;clear
%% Starting with n = 19
n = 79; % Choose any value of n from 19,39,79
h = 1/(n+1);
t = 0 : h : 1;
%% Part a
K = zeros(n+1,1);
A = zeros(n,n);
for i = 0 : n+1
    K(i+1) = 1 + (i+.5)*h;
end
%% Part b
% K = 1000;
% k_t = zeros(1,n+2);
% for m = 1 : length(t)
%     if t(m) < .5
%         k_t(m) = 1;
%     else
%         k_t(m) = K;
%     end
% end
%%
% test parameter
i = 1:n;
y = log(1 + i*h)/log(2);
b = zeros(n,1);
b(n) = K(n+1);

%% Fill up the matrix according to the given information



for j = 1 : n
    for k = 1 : n
        if j == k;
            A(j,k) = K(j) + K(j+1);
        elseif k == j-1
            A(j,k) = -K(j);
         elseif j == k-1
             A(j,k) = -K(k);           
        end
    end
end
%% LU Factorization
% [L,U] = lu_decompose(A);
 [L,U] = lu(A);
%% Now follow Algorithm from "Numerical Methhods and Computing
% Kincaid, Cheney, 6th ed, chapter 8, pg. 301:
z = zeros(n,1);
x = zeros(n,1);
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
X = x-y';
norm(X,Inf)
u = log(1+t)/log(2); %Exact solution
plot(u(1:n),'--k','LineWidth',1.2)
hold on 
% plot(x')
plot(1:n,x,'ko',1:n,x,'k-','LineWidth',1.2); %Piece-wise linear interpolation of x
set(gcf,'Color','white');
% hold(axes1,'on');

% Create ylabel
ylabel('x',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',10,...
    'Interpreter','latex');

% Create xlabel
xlabel({'n'},'Interpreter','latex','FontSize',10);
plotTickLatex2D
% Create title       
matFileName = sprintf('math609_pa1_comp_example_1_n_%d_k', n);
saveas(gcf,[matFileName,'.eps'],'epsc')
            
            







