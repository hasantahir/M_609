clf;
n = 12 :.01: 14;
lhs = 20e-7*(n.^(n+1)).*(n+1);
rhs = (2*pi).^(n+1);
plot(n,lhs,'LineWidth',1.4)
hold on
plot(n,rhs,'LineWidth',1.4)
legend('lhs','rhs')
hold off
clear
f = @(n) 1.8379*n - log((n^(n+1)*(n+1))) + 14.9603;
% f = @(n) 20e-7*n^(n+1)*(n+1) -(2*pi)^(n+1);
N = floor(fzero(f,13));
h = 2*pi/N;
cleanfigure;
matlab2tikz('filename',sprintf('transcendental1.tex'))
