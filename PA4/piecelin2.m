function v = piecelin2(x,y,u)
%Piecewise linear interpolation. (Slightly modified from Moler text)
%  v = piecelin(x,y,u) finds the piecewise linear value of S(x)
%  with S(x(j)) = y(j) and returns v = S(u).
%  First divided difference – “diff” is a built-in Matlab command
   delta = diff(y)./diff(x);
%  Find subinterval index k so that x(k) <= u < x(k+1)
   n = length(x);
   k = 1;
   for j = 2:n-1
      if x(j) <= u;  
            k = j;
      end
   end
%  Evaluate spline at u
   s = u - x(k);
   v = y(k) + s*delta(k);
