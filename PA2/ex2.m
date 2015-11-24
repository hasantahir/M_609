close all;clc;clear
N = 34; % Pick 8, 16 or 32
nd = N -2; % exclude boundaty points
n = nd^2; % number of node points
h = 1/(nd+1); % distance between points
Coeff = 4-h^2;
b = h^2*ones(n,1); %zero the excitation vector
A = zeros(n,n); % zero the matrix
TOL = 1e-12;
max_it = 100000;
w = 1.4;
options = {'Jacobi','Gauss-Seidel', 'SOR', 'SSOR'} ;
convergence = zeros(1,length(options));
norms = zeros(1,length(options));
%%
%Fill the F matrix and V vector
for i = 1:n
    for j = 1:n
        if i == j
            A(i,j) = Coeff;
        elseif j == i-1 || j == i+1
            A(i,j) = -1;
        elseif j == i-nd || j == i+nd
            A(i,j) = -1;
        end
    end
end
%
%% null out some off diagonal elements that have been set = 1 above
%
for k = nd:nd:n-nd
    A(k,k+1) = 0;
    A(k+1,k) = 0;
end
%

for iteration = 1 : length(options)
    [ x, error_norm, count] = iterative_solve(A,b,TOL,max_it,options{iteration},w);
    %
    %% Incorporate the node position vectors and potentials to include the
    %  boundary potentials in a 2-D matrix.
    
    X = zeros(N,N); %Solution matrix - In this case the first index
    %is associated with the row or Y-coordinate and the
    %second index with the column or x-coordinate.
    jt = 0;
    for i = 1:N
        for j = 1:N
            if i >1 && i<N  && j == N
                X(j,i)=0;
            elseif i>1 && i<N && j>1 && j<N
                jt=jt+1;
                X(j,i) = x(jt);
            end
        end
    end
    %% Plot and Display Results
    %
    figure(iteration)
    convergence(iteration) = count;
    hold on
    %Plot the function and the contours
    [xdg, ydg] = meshgrid(1:N, 1:N);
    surf(xdg,ydg,X)
    shading interp
    colormap parula
    % Create title
    title(['Convergence in ',num2str(count),' iterations for ', options{iteration}, ' method'],...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',10,...
        'Interpreter','latex');
    matFileName = sprintf('math609_pa2_comp_example_2_%d_n_%s', nd, options{iteration});
    saveas(gcf,[matFileName,'.eps'],'epsc')
end




