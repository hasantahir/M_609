% Numerical Mathematics and Computing, Fifth Edition
% Ward Cheney & David Kincaid
% Brooks/Cole Publ. Co.
% (c) 2003
%
% file:  rkf_ode45.m
%
% Matlab routine for Runge-Kutta-Fehlberg is called
% 'ode45'. Here we use this routine to solve the
% ODE over the interval [0,1] after setting the initial
% condition s and the error tolerance.
% The ODE is defined in the M-file 'ode45file1'.
close all
t = 0; x = 0; % Initial Conditions
errr = 1e-4; % Tolerance 
N = 1000; % Number of points used
tb = 5; % Time limit
options = odeset('RelTol',errr,'AbsTol',errr); % Set error criteria
tspan = linspace(t,tb,N);
[t,x] = ode45('ode45file2', tspan, x, options); % Solve ODE 
plot(t,(x),'LineWidth',1.4,'Color','black')
hold on
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman')
title('$i^{\prime}(t) = C E^{\prime\prime}(t) + \frac{1}{R} E^{\prime}(t) + \frac{1}{L} E(t), i(0) = 0$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');
grid on
% Create ylabel
ylabel('$i(t)$ (A)',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$t$ (seconds)',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');
% matlab2tikz('filename',sprintf('math609_pa5_example_1_2.tex'))
matlab2tikz('filename',sprintf('math609_pa5_comp_example_2_%d.tex', tb))