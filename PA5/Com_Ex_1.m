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
type ode45file2
t = 1; x = 0;
N = 100;
tb = 2;
tspan = linspace(t,tb,N);
errr = 1e-4;
options = odeset('RelTol',errr,'AbsTol',errr);
[t,x] = ode45('ode45file1', tspan, x, options);
plot(t,x,'LineWidth',1.4,'Color','black')
hold on
T = linspace(1,2,20);
Y = T.^2.*(exp(T)-exp(1));
plot(T,Y,'mo',...
                'LineWidth',2,...
                'MarkerEdgeColor','none',...
                'MarkerFaceColor',[1 0.5 1],...
                'MarkerSize',5)
legend('Numerical', 'Exact')
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman')
title('$y^{\prime}(t) = \frac{2}{t} y + t^2e^t, y(1) = 0$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');
grid on
% Create ylabel
ylabel('$y(t)$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$t$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');
% matlab2tikz('filename',sprintf('math609_pa5_example_1_2.tex'))
matlab2tikz('filename',sprintf('math609_pa5_comp_example_1_%d.tex', errr*1e4))