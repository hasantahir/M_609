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
alpha_0 = 75;
alpha_n = 72;
z = spline3_coeff_fixed(N,t,y,'fixed',alpha_0,alpha_n);
x = 0:.1:13;
S = zeros(1,length(x));
for i = 1 : length(x)
    S(i) = spline3_eval(N,t,y,z,x(i));
    if x(i) == 10
        S_10 = spline3_eval(N,t,y,z,x(i));
    end
end
close all
plot(t,y,'mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerSize',5)
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman')
grid on
hold on
plot(x,S,'LineWidth',1.4,'Color','black')
plot(10,S_10,'mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[0 0 0],...
    'MarkerSize',10)
title('Cubic Spline Interpolation of a traveling car with fixed ends',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create ylabel
ylabel('$Position (ft)$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$time (sec)$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');
legend('Data','Spline',['Position at t = 10 sec is ', num2str(S_10), ' ft'],'Location','Northwest')
matlab2tikz('filename',sprintf('math609_pa4_example_2.tex'))
