close all;clc;clear
N = 16; % Pick 8, 16 or 32
a = 2.; % 1/2 Length in the x-direction
nd = N -2; % exclude boundaty points
n = nd^2; % number of node points
h = 1/(nd+1); % distance between points
ah = 2*a/(nd-1); % distance between node points in the x-direction
alpa = 1; % alpha matrix variable
bta2 =(1.+alpa)/2.; % beta matrix variable
A = zeros(n,n); % zero the matrix
V = zeros(n,1); %zero the excitation vector
TOL = 1e-12;
max_it = 100000;
Vo = 100; % Voltage at the top of the trough
options = {'CG','SD'} ;
convergence = zeros(1,length(options));
norms = zeros(1,length(options));
%%
%Fill the A matrix and V vector
for i = 1:n
    for j = 1:n
        if i == j
            A(i,j) = 4.*bta2;
        elseif j == i-1 || j == i+1
            A(i,j) = -alpa;
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
% Fill the V vector
for k = nd:nd:n
    V(k,1) = alpa*Vo;
end

for iteration = 1 : length(options)
    [ x, error_norm, count] = iterative_solve(A,V,TOL,max_it,options{iteration});
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
    [xdg, ydg] = meshgrid(1:N, 1:N-1);
    surf(xdg,ydg,X(1:N-1,1:N))
    shading interp
    colormap jet
    axis tight
    set(gcf,'Color','white'); % Set background color to white
    set (gca,'FontName','times new roman','FontWeight','bold','FontSize',11); % Set axes fonts to
    % Create title
    title(['Convergence in ',num2str(count),' iterations for ', options{iteration}, ' method'],...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',12,...
        'Interpreter','latex');
    
    % Create ylabel
    ylabel('ny',...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',12,...
        'Interpreter','tex');
    
    % Create xlabel
    xlabel('nx',...
        'HorizontalAlignment','center',...
        'FontWeight','bold',...
        'FontSize',12,...
        'Interpreter','tex');
    matFileName = sprintf('math609_pa3_comp_example_3_%d_n_%s', N, options{iteration});
    %   MATLAB2TIKZ(FILENAME,...) or MATLAB2TIKZ('filename',FILENAME,...)
%     saveas(gcf,[matFileName,'.eps'],'epsc')
    matlab2tikz('filename',sprintf('math609_pa3_comp_example_3_%d_n_%s.tex', N, options{iteration}))
end




