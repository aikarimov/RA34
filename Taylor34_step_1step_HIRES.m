function  [x,err] =  Taylor34_step_1step_HIRES (fun,t0,h,x0)

%f0 = feval(fun,t0,x0);
% dim = length(x0);
% I = eye(dim);
dx = x0;
   dx(1) = -1.71*x0(1) + 0.43*x0(2) + 8.32*x0(3) + 0.0007;
   dx(2) = 1.71*x0(1) - 8.75*x0(2);
   dx(3) = -10.03*x0(3) + 0.43*x0(4) + 0.035*x0(5);
   dx(4) = 8.32*x0(2) + 1.71*x0(3) - 1.12*x0(4);
   dx(5) = -1.745*x0(5) + 0.43*x0(6) + 0.43*x0(7);
   dx(6) = -280*x0(6)*x0(8) + 0.69*x0(4) + 1.71*x0(5)  - 0.43*x0(6) + 0.69*x0(7);
   dx(7) = 280*x0(6)*x0(8) - 1.81*x0(7);
   dx(8) = -280*x0(6)*x0(8) + 1.81*x0(7);
  
f0 = dx;

J=[ -1.71 0.43 8.32 0 0 0 0 0; 
    1.71 -8.75 0 0 0 0 0 0;
    0 0 -10.03 0.43 0.035 0 0 0;
    0 8.32 1.71 -1.12 0 0 0 0;
    0 0 0 0 -1.745 0.43 0.43 0;
    0 0 0 0.69 1.71 (-280*x0(8)-0.43) 0.69 (-280*x0(6));
    0 0 0 0 0 (280*x0(8)) -1.81 (280*x0(6));
    0 0 0 0 0 (-280*x0(8)) 1.81 (-280*x0(6))];
Jf = J*f0;
dJf = [ 0 0 0 0 0 0 0 0; 
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 (-280*dx(8)) 0 (-280*dx(6));
    0 0 0 0 0 (280*dx(8)) 0 (280*dx(6));
    0 0 0 0 0 (-280*dx(8)) 0 (-280*dx(6))];
dJJf= [ 0 0 0 0 0 0 0 0; 
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 (-280*Jf(8)) 0 (-280*Jf(6));
    0 0 0 0 0 (280*Jf(8)) 0 (280*Jf(6));
    0 0 0 0 0 (-280*Jf(8)) 0 (-280*Jf(6))];
%j66=-280*(-280*x0(6)*x0(8) + 1.81*x0(7));
%dj66=78400*x0(8)*dx(6)+(-280*1.81*dx(7))+78400*x0(6)*dx(8);
dj66 = 78400*x0(8)*dx(6) - 506.8*dx(7) + 78400*x0(6)*dx(8);

%j68=-280*(-280*x0(6)*x0(8) + 0.69*x0(4) + 1.71*x0(5)  - 0.43*x0(6) + 0.69*x0(7));
% dj68= -280*0.69*dx(4)+(-280*1.71*dx(5))+((78400*x0(8)+0.43)*dx(6))+(-280*0.69*dx(7))+78400*x0(6)*dx(8);
%dj68= -280*0.69*dx(4)+(-280*1.71*dx(5))+((78400*x0(8)+0.43)*dx(6))+(-280*0.69*dx(7))+78400*x0(6)*dx(8);
dj68= -280*0.69*dx(4) - 280*1.71*dx(5)  - 280*(- 0.43  - 280*x0(8))*dx(6)  - 280*0.69*dx(7) + 78400*x0(6)*dx(8);


%j76=280*(-280*x0(6)*x0(8) + 1.81*x0(7)); 
%dj76=-78400*x0(8)*dx(6)+280*1.81*dx(7)-78400*x0(6)*dx(8);
%dj76=280*(-280*dx(6)*x0(8) + 1.81*dx(7) -280*x0(6)*dx(8));
dj76 = -dj66;


%j78=-280*(-280*x0(6)*x0(8) + 0.69*x0(4) + 1.71*x0(5)  - 0.43*x0(6) + 0.69*x0(7));
%dj78=280*0.69*dx(4)+(280*1.71*dx(5))+((-78400*x0(8)+0.43)*dx(6))+(280*0.69*dx(7))-78400*x0(6)*dx(8);
%dj78=280*(0.69*dx(4) + 1.71*dx(5)  + (- 0.43  - 280*x0(8))*dx(6)  + 0.69*dx(7) - 280*x0(6)*dx(8));
dj78 = -dj68;

dj86 = dj66;
dj88 = dj68;
ddJ=  [ 0 0 0 0 0 0 0 0; 
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 dj66 0 dj68;
    0 0 0 0 0 dj76 0 dj78;
    0 0 0 0 0 dj86 0 dj88];

ddf = (J^2 + dJf);
%d3f = ddJ + 2*dJ*J + J*ddf;
%d3f = ddJ + 2*dJf*J + 2*J*dJf + J^3;
d3f = ddJ + dJJf + dJf*J + 2*J*dJf + J^3;

x2 = x0 + h*f0 + 0.5*h^2*J*f0 + h^3/6*ddf*f0;%(I + 0.5*h*J + h^2/6*ddf)*x0;%3 ord
x = x2 + h^4/24*d3f*f0; %(I + 0.5*h*J + h^2/6*ddf + h^3/24*d3f)*x0;%4 ord

%************ PREV ONE ***********
% Dnum3 = Dnum4 - h^4/24*d3f*f0;
% dx = U\(L\Dnum3);
%************ ALT ONE ************


%************ ERR ESTIMATOR *************
err = x - x2;

end