% MATH 609 - Programming Assignment 1
% Computational Exercise 1
% This program solves a linear system with a Hilbert matrix
close all;clc;clear
%% Starting with n = 19
n =39; % Choose any value of n from 19,39,79
K_constant = 1000; % Choose any value from 2, 100, 1000
TOL = 1e-12;
max_it = 100000;
options = {'CG','SD'} ;
convergence = zeros(1,length(options));
norms = zeros(1,length(options));
h = 1/(n+1);
K = zeros(n+1,1);
A = zeros(n,n);
b = zeros(n,1);
t = 0 : h : 1;
%
%% Part a

% part = 'a';
% for i = 0 : n+1
%     K(i+1) = 1 + (i+.5)*h;
% end

%% Part b
% %
part = 'b';
for i = 0 : n+1
    if t(i+1) < .5
        K(i+1) = 1;
    else
        K(i+1) = K_constant;
    end
end
% %
%% Fill up the matrix A
%
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
%
%% Set up the right-hand side vector b
%
b(n) = K(n+1);
%
%% Analytical Solver for part a
%
i = 1:n;
y = log(1 + i*h)/log(2);
%% Call Iterative Solver
%
for iteration = 1 : length(options)
    [ x, error_norm, count] = iterative_solve(A,b,TOL,max_it,options{iteration});
    %
    %% Plot and Display Results
    %
    figure(iteration)
    X = x-y'; % Calculate norm
    convergence(iteration) = count;
    norms(iteration) = norm(X,2)
    u = log(1+t)/log(2); %Exact solution
    if part == 'a'
        
        plot(u(1:n),'--k','LineWidth',1.2)
        hold on
    end
    
    plot(1:n,x,'ko',1:n,x,'k-','LineWidth',1.2); %Piece-wise linear interpolation of x
    set (gca,'FontName','times new roman','FontWeight','bold','FontSize',11); % Set axes fonts to
    set(gcf,'Color','white');
    
    if part == 'a'
        
        legend({'Exact Solution', 'Numerical Solution'},'Location','northwest','FontWeight','normal')
    elseif part == 'b'
        legend({sprintf('K = %d', K_constant)},'Location','northwest','FontWeight','normal')
        
    end
    
    % Create ylabel
    ylabel('x',...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',12,...
        'Interpreter','latex');
    
    % Create xlabel
    
    xlabel('n',...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',12,...
        'Interpreter','latex');
    
    plotTickLatex2D
    % Create title
    if count < max_it
        
        title(['Convergence in ',num2str(count),' iterations for ', options{iteration}, ' method'],...
            'HorizontalAlignment','center',...
            'FontWeight','bold',...
            'FontSize',12,...
            'Interpreter','latex');
    else
        
        title(['Did not Converge in ',num2str(count),' iterations for ', options{iteration}, ' method'],...
            'HorizontalAlignment','center',...
            'FontWeight','bold',...
            'FontSize',12,...
            'Interpreter','latex');
    end
    
    matFileName = sprintf('math609_pa3_comp_example_1_n_%d_k_%s_part_%s', n, options{iteration},part);
    
    %saveas(gcf,[matFileName,'.eps'],'epsc')
    
    if part == 'a'
        
        matlab2tikz('filename',sprintf('math609_pa3_comp_example_1_n_%d_%s_part_%s.tex', n, options{iteration},part));
        
    elseif part == 'b'
        
        matlab2tikz('filename',sprintf('math609_pa3_comp_example_1_n_%d_K_%d_%s_part_%s.tex', n, K_constant, options{iteration},part));
        
    end
    
end