function [ x, error_norm, count] = iterative_solve(A,b,TOL,max_it,option)
% MATH 609, TAMU Fall 2015
% Hasan Tahir Abbas
% 10/11/2015

% This function iteratively solves the linear systerm Ax=b through
%
% 1. Steepest Descent
% 2. Conjugate-Gradient
%
%    Inputs:
%        real A(N,N): the symmetric positive definite matrix.
%        real b(N): the right hand side vector.
%        integer max_it: the maximum number of iterations.
%        real TOL: an error tolerance.
%        integer option:
% 1. SD
% 2. CG
%
%    Outputs:
%        real x(x): the solution.
%        real error_norm:  the norm of the error.
%        integer count:  the number of iterations performed.
%
%
%
%% Parameters
%
%
n = length(A);
global error_norm count x
count = 0; % Count of number of interations
x = zeros(n,1); % Initial guess of zeros
r = b - A * x; % Define residual
bnrm2 = norm( b );
error_norm = norm ( r ) / bnrm2; % Define Error norm
%
%% Definitions
%

%%
switch option
    
    case 'CG'
        for iter = 1 : max_it
            
            z = r;
            rho = ( r' * z );
            %
            %  Compute P, the direction vector.
            %
            if ( iter == 1)
                p = z;
            else
                beta = rho / rho_1;
                p = z + beta * p;
            end
            
            q = A * p;
            alpha = rho / ( p' * q );
            %
            %  Update the approximate solution.
            %
            x = x + alpha * p;
            %
            %  Compute the residual.
            %
            r = r - alpha * q;
            error_norm = norm ( r ) / bnrm2;
            count = iter;
            
            if ( error_norm <= TOL )
                break
            end
            
            rho_1 = rho;
            
        end
        
    case 'SD'
        
        
        
        for iter = 1 : max_it
            
            v = b - A * x;
            q = A*v;
            t = (v' * v) / (v' * q );
            x = x + t*v;
            error_norm = norm(v)/bnrm2;
            if ( error_norm <= TOL )
                break
            end
            
        end
        count = iter;
       
        
end

end
