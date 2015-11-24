function z = spline3_coeff(N,t,y,options,alpha_0,alpha_n)
% function z = spline3_coeff(n,t,y,data);
% Based on the pseudocode spline3_coeff presented in Kincaid and Cheney,
% " Numerical Mathematics and Computing"
%
% The natural spline for interpolating data at the knots a=x0<...<x_n=b;
% System of equations is solved using the tridiagonal nature
%    WARNING: THIS PROGRAM REQUIRES AT LEAST 4 KNOTS
%  INPUT:
%    N - number of knots
%    t  - knot points
%    y  - interpolation points
%    options - free or fixed ends
%  OUTPUT:
%     z - column vector containing solution of the triadiagonal matrix
%         equaion
                                                                               
%% determine lengths of the intervals
n = N - 1;
h = zeros(1,n);
b = zeros(1,n);
u = zeros(1,n);
v = zeros(1,n);
z = zeros(1,n+1);
for i = 1 : n
   
    h(i) = t(i+1) - t(i);
    b(i) = (y(i+1) - y(i))/h(i);

end

u(1) = 2*(h(1) + h(2));  % 2(h(0) + h(1))
v(1) = 6*(b(2) - b(1));  % 6(b(1) - b(0))
%% determine the right hand side of the system (without first and last eqn)
for i = 2 : n
    
    u(i) = 2*(h(i) + h(i-1)) - h(i-1)^2/u(i-1);  % 2(h(0) + h(1))
    v(i) = 6*(b(i) - b(i-1)) - h(i-1)*v(i-1)/u(i-1);  % 6(b(1) - b(0))
    
end   
%% Solve the Tridiagonal System
switch options
    
    case 'free'
        z(1) = 0;
        z(n+1) = 0;
        for i = n : -1 : 2
            
            z(i) = (v(i) - h(i)*z(i+1))/u(i);
            
        end
        
    case 'fixed'
        for i = n : -1 : 2
            
            z(i) = (v(i) - h(i)*z(i+1))/u(i);
            
        end
        H = t(2) - t(1);
        
        z(1) = (alpha_0 + z(2)*H/6 + y(1)/H - y(2)/H)*6/(4*H);
        
        H = t(n) - t(n-1);
        z(n+1) = (alpha_n - z(n)*H/6 - (y(n)- y(n-1))/H)*3/H;
        
    otherwise
        disp('Error: Wrong Option inserted');
end
end