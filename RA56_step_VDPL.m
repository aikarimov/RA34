function  [x,err, J2, denCond] =  RA56_step_VDPL (fun,t0,h,x0,mu)

dx = x0;

x12 = x0(1)^2;
mu2 = 2*mu;
x1x2 = x0(1)*x0(2);

dx(1) = x0(2);
dx(2) = mu*(1 - x12)*x0(2) - x0(1);

F = [dx(1);dx(2)];

f0 = feval(fun,t0,x0);
dim = length(x0);
I = eye(dim);

j21 = - mu2*x1x2 - 1;
j22 = mu*(1 - x12);
J = [ 0 1; j21 j22];

JF = J*f0;

dj21x= - mu2*x0(2);
dj21y= - mu2*x0(1);
dj22x= - mu2*x0(1);
dj22y= 0;
dJF  = [ 0, 0; [dj21x,dj21y]*dx, [dj22x,dj22y]*dx ];
dJJF = [ 0, 0; [dj21x,dj21y]*JF, [dj22x,dj22y]*JF ];

dJFF = dJF*F;

dJdJFF = [ 0, 0; [dj21x,dj21y]*dJFF, [dj22x,dj22y]*dJFF ];

JJF = J*JF;
JJJF = J*JJF;

dJJJF =  [ 0, 0; [dj21x,dj21y]*JJF, [dj22x,dj22y]*JJF ];
dJJJJF = [ 0, 0; [dj21x,dj21y]*JJJF, [dj22x,dj22y]*JJJF ];

dJdJJFF = [ 0, 0; [dj21x,dj21y]*(dJJF*F), [dj22x,dj22y]*(dJJF*F) ];
dJdJFJF = [ 0, 0; [dj21x,dj21y]*(dJF*JF), [dj22x,dj22y]*(dJF*JF) ];
dJJdJFF = [ 0, 0; [dj21x,dj21y]*(J*dJF*F), [dj22x,dj22y]*(J*dJF*F) ];

ddJFF =  [0, 0;-4*mu*dx(1)*dx(2) - mu2*dx(1)*dx(2), -2*mu*dx(1)^2];

dJddJFFF = [ 0, 0; [dj21x,dj21y]*(ddJFF*F), [dj22x,dj22y]*(ddJFF*F)];

ddJJFF = [0, 0;-4*mu*JF(1)*dx(2) - mu2*JF(1)*dx(2), -2*mu*JF(1)*dx(1)];
ddJFJF = [0, 0;-4*mu*dx(1)*JF(2) - mu2*dx(1)*JF(2), -2*mu*dx(1)*JF(1)];
ddJdJFFF = [0, 0;-4*mu*dJFF(1)*dx(2) - mu2*dJFF(1)*dx(2), -2*mu*dJFF(1)*dx(1)];
ddJJJFF =  [0, 0;-4*mu*JJF(1)*dx(2) - mu2*JJF(1)*dx(2), -2*mu*JJF(1)*dx(1)];
ddJJFJF =  [0, 0;-4*mu*JF(1)*JF(2) - mu2*JF(1)*JF(2), -2*mu*JF(1)*JF(1)];
ddJFdJFF = [0, 0;-4*mu*dx(1)*dJFF(2) - mu2*dx(1)*dJFF(2), -2*mu*dx(1)*dJFF(1)];
ddJFJJF  = [0, 0;-4*mu*dx(1)*JJF(2) - mu2*dx(1)*JJF(2), -2*mu*dx(1)*JJF(1)];

ddf = (J^2 + dJF);
d3f = ddJFF + dJJF + 2*dJF*J + J*dJF + J^3;
d4f = ddJJFF + 2*ddJFJF+ 3*ddJFF*J + dJdJFF + dJJJF + 3*dJJF*J + 3*dJF*dJF + 3*dJF*J*J + J*ddJFF + J*dJJF + 2*J*dJF*J + J*J*dJF + J^4;
d5f = ddJdJFFF + ddJJJFF + 3*ddJJFJF + 4*ddJJFF*J + 3*ddJFdJFF + 3*ddJFJJF + 8*ddJFJF*J + 6*ddJFF*dJF + 6*ddJFF*J*J + dJddJFFF + dJdJJFF + 2*dJdJFJF + 4*dJdJFF*J + dJJdJFF + dJJJJF + 4*dJJJF*J + 6*dJJF*dJF + 6*dJJF*J*J + 4*dJF*ddJFF + 4*dJF*dJJF + 8*dJF*dJF*J + 4*dJF*J*dJF + 4*dJF*J*J*J + J*ddJJFF + 2*J*ddJFJF + 3*J*ddJFF*J + J*dJdJFF + J*dJJJF + 3*J*dJJF*J + 3*J*dJF*dJF + 3*J*dJF*J*J + J*J*ddJFF + J*J*dJJF + 2*J*J*dJF*J + J*J*J*dJF + J^5;

%+ h^4/360*( 6*d4f + 10*ddf^2 - 15*d3f*J)
Dnum6 = (I + h^2*(-1/4*J^2 + 1/3*ddf) + h^4/360*( 6*d4f + 10*ddf^2 - 15*d3f*J) )*h*f0;
Den6 = I - 0.5*h*J + h^2/6*ddf - h^3/24*d3f + h^4/120*d4f - h^5/720*d5f;

denCond = cond(Den6);

[L,U] = lu(Den6);
dx = U\(L\Dnum6);
x = x0 + dx; %6 ord

%Dnum5 = Dnum6 - h^6/720*d5f*f0;
%dx = U\(L\Dnum5);
Den5 = I - 0.5*h*J + h^2/6*ddf - h^3/24*d3f + h^4/120*d4f;
dx = Den5\Dnum6;

x2 = x0 + dx; %5 ord

J2=J;
err = x - x2;
end
