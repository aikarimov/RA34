function  [x,err] =  Taylor34_step_1step_RICCATI (fun,t0,h,x0,aparam)

dx = x0;
dx(1) = aparam - x0(1)^2 - x0(3)*x0(2);
dx(2) = - x0(2)*(x0(1) + x0(4));
dx(3) = - x0(3)*(x0(1) + x0(4));
dx(4) = aparam - x0(4)^2 - x0(3)*x0(2);
  
f0 = dx;

j11=-2*x0(1);
j12=-x0(3);
j13=-x0(2);
%j21=j13;
j22=-x0(1)-x0(4);
%j24=j21
%j31=j12;
%j33=j22;
%j34=j31
% j42=j12;
% j43=j13;
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

dj11= 4*x0(1)*dx(1)+2*x0(3)*dx(2)+2*x0(2)*dx(3);

dj12= x0(3)*dx(1)+(x0(1) + x0(4))*dx(3)+x0(3)*dx(4);

dj13= x0(2)*dx(1)+(x0(1) + x0(4))*dx(2)+x0(2)*dx(4);

dj22= 2*x0(1)*dx(1)+2*x0(3)*dx(2)+2*x0(2)*dx(3)+2*x0(4)*dx(4);
ddJff =  [dj11 dj12 dj13 0;
       dj13 dj22 0 dj13
       dj12 0 dj22 dj12
       0 dj12 dj13 dj11];

ddf = (J^2 + dJf);

d3f = ddJff + dJJf + 2*dJf*J + J*dJf +  J^3;

x2 = x0 + h*f0 + 0.5*h^2*J*f0 + h^3/6*ddf*f0;
x = x2 + h^4/24*d3f*f0;

err = x - x2;

end