function [x, err, denCond] = RA34elwise_step_VDPL(fun,t0,h,x0,mu)


dim = length(x0);
I = eye(dim);
f0 = feval(fun,t0,x0);

dx = x0;
dx(1) = x0(2);
dx(2) = mu*(1 - x0(1)^2)*x0(2) - x0(1);

j21=-2*mu*x0(1)*x0(2)-1;
j22=mu-mu*x0(1)*x0(1);
J=[ 0 1; j21 j22];

dj21x=-2*mu*x0(2);
dj21y=-2*mu*x0(1);
dj22x=-mu*2*x0(1);
dj22y= 0;
dJ= [ 0 0; (dj21x*dx(1)+dj21y*dx(2)) (dj22x*dx(1)+dj22y*dx(2)) ];

ddj21x=0+(-2*mu*mu*x0(2)+3*mu*x0(1)*x0(1)*x0(2)+4*mu*x0(1));
ddj21y=-2*mu*2*x0(2)-2*mu*x0(1)*mu + 2*mu*x0(1)*mu*x0(1)^2;
ddj22x=-mu*2*x0(2);
ddj22y= -mu*2*x0(1);
ddJ= [ 0 0; (ddj21x*dx(1)+ddj21y*dx(2)) (ddj22x*dx(1)+ddj22y*dx(2)) ];

ddf = (J^2 + dJ);
d3f = ddJ + 2*dJ*J + J*ddf;
Dnum4 = ((I + h^2*(-1/4*J^2 + 1/3*ddf))*h*f0).*f0;
Den4 = (I - 0.5*h*J + h^2/6*ddf - h^3/24*d3f);

denCond = cond(Den4);
Den4 = Den4*f0;


dx = Dnum4./Den4;
x = x0 + dx; %4 ord

Dnum3 = Dnum4 - h^4/24*d3f*f0.*f0;
dx = Dnum3./Den4;
x2 = x0 + dx; %3 ord
err = x - x2;

end