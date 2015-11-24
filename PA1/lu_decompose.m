% This program computes the LU Decomposition of a matrix A of size n
% It also assumes the factorization exists without restructuring the
% matrix

% It is based on the row-based Doolittle Algorithm (all l_kk = 1)

% Hasan Tahir Abbas
% 9/10/2015

function [L,U] = lu_decompose(A)
% where,
% L is the lower triangular matrix
% U is the upper triangular matrix
% A is the matrix to be factorized
%

n = length(A); % Size of the matrix
L = eye(n);
U = zeros(n);
temp_sum_u = 0;
temp_sum_l = 0;

for k = 1 : n
    for j = k+1 : n
        for s = 1 : k-1
            temp_sum_u = temp_sum_u + L(k,s)*U(s,j);
        end
        U(k,j) = A(k,j) - temp_sum_u;
    end
    
    for i = k+1 : n
        for s = 1 : k-1
            temp_sum_l = temp_sum_l + L(i,s)*U(s,k);
        end
        L(k,j) = (A(k,j) - temp_sum_u)/U(k,k);
    end
end
