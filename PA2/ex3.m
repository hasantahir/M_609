clc; close all;clear %clears the command window and closes the plot window when the code is rerun
% This is a program to solve for the potential in a square cylinder
%%
% Essential coefficients and constants
Vo = 100.; % Voltage potential
a = 2.; % 1/2 Length in the x-direction
b = 3.; % Total length in the y-direction
nxd = 10; % number of data points in the x-direction
nyd = 15; % number of data points in the y-direction
nx=nxd-2; % computed data points - without boundary points
ny=nyd-2;
n = nx*ny; % number of node points
ah = 2*a/(nxd-1); % distance between node points in the x-direction
bk = b/(nyd-1); % distance between node points in the y-direction
alpa = ah/bk; % alpha matrix variable
alpa2 = alpa^2;
bta2 =(1.+alpa2)/2.; % beta matrix variable
V = zeros(n,1); %zero the excitation vector
F = zeros(n,n); % zero the matriz
TOL = 1e-12;
max_it = 100000;
w = 1.4;
options = {'Jacobi','Gauss-Seidel', 'SOR', 'SSOR'} ;
convergence = zeros(1,length(options))';
norms = zeros(1,length(options));
%%
%Fill the F matrix and V vector
for i = 1:n
    for j = 1:n
        if i == j
            F(i,j) = 4.*bta2;
        elseif j == i-1 || j == i+1
            F(i,j) = -alpa2;
        elseif j == i-ny || j == i+ny
            F(i,j) = -1.;
        end
    end
end
% null out some off diagonal elements that have been set = alpa2 above
for k = ny:ny:n-ny
    F(k,k+1) = 0.;
    F(k+1,k) = 0.;
end
% Fill the V vector
for k = ny:ny:n
    V(k,1) = alpa2*Vo;
end
for iteration = 1 : length(options)
    [ x, error_norm, count] = iterative_solve(F,V,TOL,max_it,options{iteration},w);
    %
    %% Incorporate the node position vectors and potentials to include the
    %  boundary potentials in a 2-D matrix.
    
    X = zeros(nyd,nxd); %Solution matrix - In this case the first index
    %is associated with the row or Y-coordinate and the
    %second index with the column or x-coordinate.
    jt = 0;
    for i = 1:nxd
        for j = 1:nyd
            if i >1 && i<n  && j == n
                X(j,i)=0;
            elseif i>1 && i<nxd && j>1 && j<nyd
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
    yd = 0:bk:b;
    xd = -a:ah:a;
    [xdg, ydg] = meshgrid(xd, yd);
    surf(xdg,ydg,X)
    shading interp
    colormap parula
    % Create title
    title(['Convergence in ',num2str(count),' iterations for ', options{iteration}, ' method'],...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',10,...
        'Interpreter','latex');
    matFileName = sprintf('math609_pa3_comp_example_2_%d_n_%s', n, options{iteration});
    saveas(gcf,[matFileName,'.eps'],'epsc')
end





