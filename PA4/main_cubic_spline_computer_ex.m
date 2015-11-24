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
filename = 'computer_ex.txt';
file_read = importdata(filename);
x = file_read(:,1); % the first column has the t vector of length N+1
y = file_read(:,2); % the second column has the y vector of length N+1
N = length(x); % Number of knots
t = 1 : N;
len = 101;
z1 = spline3_coeff(N,t,x,'free');
x1 = linspace(min(t),max(t),len);
S1 = zeros(1,length(x1));
for i = 1 : length(x1)
    S1(i) = spline3_eval(N,t,x,z1,x1(i));
end


z2 = spline3_coeff(N,t,y,'free');
x2 = linspace(min(t),max(t),len);
S2 = zeros(1,length(x2));
for i = 1 : length(x2)
    S2(i) = spline3_eval(N,t,y,z2,x2(i));
end
close all
%% Plotting Routine
%
% Plot of the x-spline
figure
plot(t,x,'mo',...
                'LineWidth',2,...
                'MarkerEdgeColor','none',...
                'MarkerFaceColor',[0.5 0.5 0.5],...
                'MarkerSize',5)
hold on
plot(x1,S1,'LineWidth',1.4,'Color','black')
legend('Data','Spline')
%
% Plot of the y-spline
figure
plot(t,y,'mo',...
                'LineWidth',2,...
                'MarkerEdgeColor','none',...
                'MarkerFaceColor',[0.5 0.5 0.5],...
                'MarkerSize',5)
hold on
plot(x2,S2,'LineWidth',1.4,'Color','black')
legend('Data','Spline')
%
% Plot of the space-spline
figure
plot(x,y,'mo',...
                'LineWidth',2,...
                'MarkerEdgeColor','none',...
                'MarkerFaceColor',[0.5 0.5 0.5],...
                'MarkerSize',5)
hold on
plot(S1,S2,'LineWidth',1.4,'Color','black')
legend('Data','Spline')
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman')
grid on
title(['Parametric Cubic Spline for a hand-drawn curve with n = ',num2str(N)],...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create ylabel
ylabel('$y$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$x$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');
matlab2tikz('filename',sprintf('math609_pa4_example_3_2.tex'))