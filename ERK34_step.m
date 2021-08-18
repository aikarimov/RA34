function [x, err] = ERK34_step(fun,t0,h,x0)
%
b1 = 1/6;
b2 = 1/3;
b3 = 1/3;
b4 = 1/6;

lam = 0.1;
bz1 = 1/6;
bz2 = 1/3;
bz3 = 1/3;
bz4 = 1/6 - lam;
bz5 = lam;

k1 = feval(fun,t0,x0);
x2 = x0 + 0.5*h*k1;
k2 = feval(fun,t0 + 0.5*h, x2);
x3 = x0 + 0.5*h*k2;
k3 = feval(fun,t0 + 0.5*h, x3);
x4 = x0 + h*k3;
k4 = feval(fun,t0 + h, x4);
x5 = x0 + h*(b1*k1 + b2*k2 + b3*k3 + b4*k4);
k5 = feval(fun,t0 + h, x5);

x  = x5; %4 ord
x2 = x0 + h*(bz1*k1 + bz2*k2 + bz3*k3 + bz4*k4 + bz5*k5); %3 ord

err = x - x2;

end