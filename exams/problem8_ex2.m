clc

syms A1 A2 x1 x2
S = solve([A1 + A2 == 1/2,...
    A1*x1 + A2*x2 == 1/3,...
    A1*x1^2 + A2*x2^2 == 1/4,...
    A1*x1^3 + A2*x2^3 == 1/5 ],...
    [A1, A2, x1, x2]);
a1 = S.A1
a2 = S.A2
x1 = S.x1
x2 = S.x2

eval(a1)
eval(a2)
eval(x1)
eval(x2)


eval(a1(1)*sin(x1(1)) + a2(1)*sin(x2(1)))
eval(a1(2)*sin(x1(2)) + a2(2)*sin(x2(2)))

eval(a1(1)*cos(x1(1)) + a2(1)*cos(x2(1)))
eval(a1(2)*cos(x1(2)) + a2(2)*cos(x2(2)))
