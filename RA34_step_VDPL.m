function  [x,err, J2, denCond] =  RA34_step_VDPL (fun,t0,h,x0,mu)

dx = x0;

x12 = x0(1)^2;
mu2 = 2*mu;
x1x2 = x0(1)*x0(2);

dx(1) = x0(2);
dx(2) = mu*(1 - x12)*x0(2) - x0(1);
f0 = feval(fun,t0,x0);
dim = length(x0);
I = eye(dim);

j21 = - mu2*x1x2 - 1;
j22 = mu*(1 - x12);
J = [ 0 1; j21 j22];

Jf = J*f0;

dj21x= - mu2*x0(2);
dj21y= - mu2*x0(1);
dj22x= - mu2*x0(1);
dj22y= 0;
dJf= [ 0 0; (dj21x*dx(1)+dj21y*dx(2)) (dj22x*dx(1)+dj22y*dx(2)) ];
dJJf = [ 0 0; (dj21x*Jf(1)+dj21y*Jf(2)) (dj22x*Jf(1)+dj22y*Jf(2)) ];

ddJff = [0, 0;-4*mu*dx(1)*dx(2) - mu2*dx(1)*dx(2), -2*mu*dx(1)^2];

ddf = (J^2 + dJf);
d3f = ddJff + dJJf + 2*dJf*J + J*dJf + J^3;


Dnum4 = (I + h^2*(-1/4*J^2 + 1/3*ddf))*h*f0;
Den4 = I - 0.5*h*J + h^2/6*ddf - h^3/24*d3f;

denCond = cond(Den4);

[L,U] = lu(Den4);
dx = U\(L\Dnum4);
x = x0 + dx; %4 ord

Dnum3 = Dnum4 - h^4/24*d3f*f0;
dx = U\(L\Dnum3);

% Den3 = I - 0.5*h*J + h^2/6*ddf;
% dx = Den3\Dnum4;


x2 = x0 + dx; %3 ord

J2=J;
err = x - x2;
end