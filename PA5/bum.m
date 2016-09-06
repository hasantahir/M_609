clear all
h = 0.0001;       % step size
t = -0:h:5;    % domain
E = exp(-.06*pi*t).*sin(2*t-pi);      % range
E1 = gradient(E)/h;   % first derivative
E2 = gradient(E1)/h;   % second derivative
plot(t,E1,'r',t,E,'b', t,E2,'k')