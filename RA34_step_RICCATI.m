function  [x,err, J2, denCond] =  RA34_step_RICCATI (fun,t0,h,x0,a)

dim = length(x0);
I = eye(dim);

dx = x0;
dx(1) = a - x0(1)^2 - x0(3)*x0(2);
dx(2) = - x0(2)*(x0(1) + x0(4));
dx(3) = - x0(3)*(x0(1) + x0(4));
dx(4) = a - x0(4)^2 - x0(3)*x0(2);
  
f0 = dx;

j11=-2*x0(1);
j12=-x0(3);
j13=-x0(2);
j22=-x0(1)-x0(4);

j44=-2*x0(4);

J=[j11 j12 j13 0;
   j13 j22 0 j13;
   j12 0 j22 j12;
   0 j12 j13 j44];

Jf = J*f0;
dJf = [-2*dx(1) -1*dx(3) -1*dx(2) 0;
      -1*dx(2) (-1*dx(1)-1*dx(4)) 0 -1*dx(2);
      -1*dx(3) 0 (-1*dx(1)-1*dx(4)) -1*dx(3);
      0 -1*dx(3) -1*dx(2) -2*dx(1);];
dJJf = [-2*Jf(1) -1*Jf(3) -1*Jf(2) 0;
      -1*Jf(2) (-1*Jf(1)-1*Jf(4)) 0 -1*Jf(2);
      -1*Jf(3) 0 (-1*Jf(1)-1*Jf(4)) -1*Jf(3);
      0 -1*Jf(3) -1*Jf(2) -2*Jf(1);];  

% dj11= 4*x0(1)*dx(1)+2*x0(3)*dx(2)+2*x0(2)*dx(3);
% 
% dj12= x0(3)*dx(1)+(x0(1) + x0(4))*dx(3)+x0(3)*dx(4);
% 
% dj13= x0(2)*dx(1)+(x0(1) + x0(4))*dx(2)+x0(2)*dx(4);
% 
% dj22= 2*x0(1)*dx(1)+2*x0(3)*dx(2)+2*x0(2)*dx(3)+2*x0(4)*dx(4);
% ddJff=  [dj11 dj12 dj13 0;
%        dj13 dj22 0 dj13
%        dj12 0 dj22 dj12
%        0 dj12 dj13 dj11];

ddJff = zeros(4);

ddf = (J^2 + dJf);

d3f = ddJff + dJJf + 2*dJf*J + J*dJf +  J^3;

Dnum4 = (I + h^2*(-1/4*J^2 + 1/3*ddf))*h*f0;
Den4 = I - 0.5*h*J + h^2/6*ddf - h^3/24*d3f;

denCond = cond(Den4);

[L,U] = lu(Den4);
dx = U\(L\Dnum4); 
dx = Den4\Dnum4;
x = x0 + dx; %4 ord

Dnum3 = Dnum4 - h^4/24*d3f*f0;
dx = U\(L\Dnum3);
x2 = x0 + dx; %3 ord
J2=J;
err = x - x2;

end