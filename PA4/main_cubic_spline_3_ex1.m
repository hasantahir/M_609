% Main program to evaluate cubic spline
% Hasan Tahir Abbas
% 01/11/2015
% Math 609 : Programming Assignment 4
%
%
% This program reads data files and interpolates using cubic spline
% interpolation
clear;close all
%% Data
n = 21;
t = linspace(0,pi,n+1);
T = linspace(0,pi,101);
x = cos(t);
x1 = -sin(t);
x2 = -cos(t);
y = sin(t);
y1 = cos(t);
y2 = -sin(t);
%% cos(t)
z1 = spline3_coeff(n+1,t,x,'free');
S_1 = zeros(1,length(T));
S1_1 = zeros(1,length(T));
S2_1 = zeros(1,length(T));
for i = 1 : length(T)
    [S_1(i),S1_1(i),S2_1(i)] = spline3_eval_multi(n,t,x,z1,T(i));
end
%% Sin(t)
z2 = spline3_coeff(n+1,t,y,'free');
S_2 = zeros(1,length(T));
S1_2 = zeros(1,length(T));
S2_2 = zeros(1,length(T));
for i = 1 : length(T)
    [S_2(i),S1_2(i),S2_2(i)] = spline3_eval_multi(n,t,y,z2,T(i));
end
close all
plot(x,y,'mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerSize',5)
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman')
hold on
plot(S_1,S_2,'LineWidth',1.4,'Color','black')
axis([-2 2 0 1.5])
grid on
title(['Parametric Cubic Spline for x(t) = cos(t), y(t) = sin(t) with $0 <= t <=\pi$ and n = ',num2str(n)],...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create ylabel
ylabel('$sin(t)$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$cos(t)$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');
legend('Data','Spline')
matlab2tikz('filename',sprintf('math609_pa4_example_3_1__n_%d.tex', n))
%% Calculation of Errors
x = cos(T);
x1 = -sin(T);
x2 = -cos(T);
y = sin(T);
y1 = cos(T);
y2 = -sin(T);
E = max(abs(x-S_1) + abs(y-S_2))
E1 = max(abs(x1-S1_1) + abs(y1-S1_2))
E2 = max(abs(x2-S2_1) + abs(y2-S2_2))