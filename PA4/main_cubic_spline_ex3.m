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
filename = 'US_population.txt';
file_read = importdata(filename);
t = file_read(:,1); % the first column has the t vector of length N+1
y = file_read(:,2); % the second column has the y vector of length N+1
N = length(t); % Number of knots
z = spline3_coeff(N,t,y,'free');
% x = [1965,1975,1985];
x = 1930: .5 : 1990;
S = zeros(1,length(x));
for i = 1 : length(x)
    S(i) = spline3_eval(N,t,y,z,x(i));
    if x(i) == 1965
        S_65 = spline3_eval(N,t,y,z,x(i));
    elseif x(i) == 1975
        S_75 = spline3_eval(N,t,y,z,x(i));
    elseif x(i) == 1985
        S_85 = spline3_eval(N,t,y,z,x(i));
    end
end
close all
plot(t,y,'mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerSize',5);
hold on
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman')
grid on
plot(x,S,'LineWidth',1.4,'Color','black')
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
plot(1965,S_65,'ms',...
    'LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[0 0 0],...
    'MarkerSize',10)

plot(1975,S_75,'md',...
    'LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[0 0 0],...
    'MarkerSize',10)

plot(1985,S_85,'mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','none',...
    'MarkerFaceColor',[0 0 0],...
    'MarkerSize',10)


% title('Cubic Spline Interpolation of US Population',...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'Interpreter','latex');

% Create ylabel
ylabel('Population (in K)',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('Year',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');
legend('Data','Spline',...
    ['Population in 1965 is ' num2str(ceil(S_65))],...
    ['Population in 1975 is ' num2str(ceil(S_75))],...
    ['Population in 1985 is ' num2str(ceil(S_85))],...
    'Location','Northwest')
matlab2tikz('filename',sprintf('math609_pa4_example_3.tex'))