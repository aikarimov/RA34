function  [x,err] =  Taylor34_step_1step_VDPL (fun,t0,h,x0,mu)

dx = x0;
dx(1) = x0(2);

x12 = x0(1)^2;
mu2 = 2*mu;
x1x2 = x0(1)*x0(2);


dx(1) = x0(2);
dx(2) = mu*(1 - x12)*x0(2) - x0(1);
f0 = feval(fun,t0,x0);


j21= - mu2*x1x2 - 1;
j22= mu*(1 - x12);
J=[ 0 1; j21 j22];

Jf = J*f0;

dj21x= - mu2*x0(2);
dj21y= - mu2*x0(1);
dj22x= - mu2*x0(1);
dj22y= 0;
dJf= [ 0 0; (dj21x*dx(1)+dj21y*dx(2)) (dj22x*dx(1)+dj22y*dx(2)) ];

dJJf = [ 0 0; (dj21x*Jf(1)+dj21y*Jf(2)) (dj22x*Jf(1)+dj22y*Jf(2)) ];


% ddj21x=0+(- mu2*mu*x0(2) + 3*mu*x12*x0(2) + 4*mu*x0(1));
% ddj21y= -4*mu*x0(2) - mu2*x0(1)*mu + mu2*x0(1)*mu*x12;
% ddj22x= - mu2*x0(2);
% ddj22y= - mu2*x0(1);
% ddJff= [ 0 0; (ddj21x*dx(1)+ddj21y*dx(2)) (ddj22x*dx(1)+ddj22y*dx(2)) ];

ddJff = [0, 0;-4*mu*dx(1)*dx(2) - mu2*dx(1)*dx(2), -2*mu*dx(1)^2];

ddf = (J^2 + dJf);
d3f = ddJff + dJJf + 2*dJf*J + J*dJf + J^3;

x2 = x0 + h*f0 + 0.5*h^2*J*f0 + h^3/6*ddf*f0;
x = x2 + h^4/24*d3f*f0;

%************ ERR ESTIMATOR *************
err = x - x2;

end