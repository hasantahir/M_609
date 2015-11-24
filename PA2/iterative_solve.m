function [ x, error_norm, count] = iterative_solve(A,b,TOL,max_it,option,varargin)
% MATH 609, TAMU Fall 2015
% Hasan Tahir Abbas
% 09/28/2015

% This function iteratively solves the linear systerm Ax=b through
%
% 1. Jacobi
% 2. Gauss-Siedel
% 3. SOR 
% 4. SSOR
%
%    Inputs:
%        real A(N,N): the symmetric positive definite matrix.
%        real b(N): the right hand side vector.
%        integer max_it: the maximum number of iterations.
%        real TOL: an error tolerance.
%        integer option: 
                % 1. Jacobi
                % 2. Gauss-Siedel
                % 3. SOR 
                % 4. SSOR
%
%    Outputs:
%        real x(x): the solution.
%        real error_norm:  the norm of the error.
%        integer count:  the number of iterations performed.
%
% The Splitting of A is done by
%
% A = D - L - U
%
%
%% Splitting of matrix A
%
L = -tril(A,-1);
D = diag(diag(A));
U = -triu(A,1);
%
%
%% Parameters
%
%
w = varargin{1};
n = length(A);
global error_norm count x
error_norm = 0;
count = 0;
x = zeros(n,1);
%
%%
switch option

    case 'Jacobi'
        %% JACOBI ITERATION
        %
        % Q = D
        % x(k) = inv(Q)*(L+U)*x(k-1) + inv(Q)*b Iteration Method
        %
        %
        % Initial Guess
        %
        Jacobi_sol = zeros(n,1); % Zeros
        Q = D;
        %
        % Start Iteration
        for iter = 1 : max_it
            
            x_1 = Jacobi_sol;
            Jacobi_sol = Q \ ((L+U)*Jacobi_sol + b);  %  Update the approximation.
            error_norm = norm(Jacobi_sol - x_1)/norm(Jacobi_sol); %Compute the error.
            
            if ( error_norm <= TOL ) %Check for convergence
                break
            end
        end
        count = iter;
        x = Jacobi_sol;
        

    case 'Gauss-Seidel'
        %% GAUSS-SEIDEL ITERATION
        %
        % Q = (D-L)^(-1)
        % x(k) = inv(Q)*(L+U)*x(k-1) + inv(Q)*b Iteration Method
        %
        %
        % Initial Guess
        %
        GS_sol = zeros(n,1); % Zeros
        Q = D - L;
        %
        % Start Iteration
        for iter = 1 : max_it
            
            x_1 = GS_sol;
            GS_sol = Q \ (U*GS_sol + b);  %  Update the approximation.
            error_norm = norm(GS_sol - x_1)/norm(GS_sol); %Compute the error.
            
            if ( error_norm <= TOL ) %Check for convergence
                break
            end
        end
        count = iter;
        x = GS_sol;
        
    case 'SOR'            
        %% SOR ITERATION
        %
        % Q = (D-w*L)
        % x(k) = inv(Q)*(w*U+(1-w)*D)*x(k-1) + inv(Q)*w*b Iteration Method
        %
        %
        % Initial Guess
        %
        SOR_sol = zeros(n,1); % Zeros
        Q = D - w*L;
        %
        % Start Iteration
        for iter = 1 : max_it
            
            x_1 = SOR_sol;
            SOR_sol = Q \ ((w*U+(1-w)*D)*SOR_sol + w*b);  %  Update the approximation.
            error_norm = norm(SOR_sol - x_1)/norm(SOR_sol); %Compute the error.
            
            if ( error_norm <= TOL ) %Check for convergence
                break
            end
        end
        count = iter;
        x = SOR_sol;

    case 'SSOR'       
        %% SSOR ITERATION
        %
        % B1 = (D-w*U)^(-1)*(w*L + (1-w)*D)
        % B2 = (D-w*L)^(-1)*(w*U + (1-w)*D)
        % x(k) = B1*B2*x(k-1) + w*(2-w)*(D-w*U)^(-1)*D*(D-w*L)^(-1)*b Iteration Method
        %
        %
        % Initial Guess
        %
        SSOR_sol = zeros(n,1); % Zeros
        DU = D-w*U;
        DL = D-w*L;
        B1 = DU\(w*L + (1-w)*D);
        B2 = DL\(w*U + (1-w)*D);
        %
        % Start Iteration
        for iter = 1 : max_it
            
            x_1 = SSOR_sol;
            SSOR_sol = B1*B2*SSOR_sol + w*(2-w)*DU\D*DL\b;  %  Update the approximation.
            error_norm = norm(SSOR_sol - x_1)/norm(SSOR_sol); %Compute the error.
            
            if ( error_norm <= TOL ) %Check for convergence
                break
            end
        end
        count = iter;
        x = SSOR_sol;
        
    otherwise 
        %% Else show error
        error(' Please enter correct Iteration option');
end
x;
error_norm;
count;
end
