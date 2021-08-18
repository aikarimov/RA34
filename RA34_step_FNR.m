function  [x,err, J2, denCond] =  RA34_step_FNR(fun,t0,h,x0)
%f0 = feval(fun,t0,x0);
dim = length(x0);
I = eye(dim);

a=[0.5,0.7,0.8,0.002,-1]';
A=0.083; w = 0.00003*2*3.14;

dx = x0;
i= x0(5);
dx(1)=x0(1)-(x0(1)*x0(1)*x0(1))/3-x0(2)-x0(3)+i;
dx(2)=a(1)*(a(2)+x0(1)-a(3)*x0(2));
dx(3)=a(4)*(a(5)-x0(3)-x0(1));
dx(4)=A*x0(5);
dx(5)=-w^2*x0(4);
 
f0 = dx;

j11= 1-x0(1)*x0(1);
j12=-1;
j15=1;
j21=a(1);
j22=-a(3)*a(1);
j31=-a(4);
j45=A;
j54=-w^2;
J=[j11 j12 j12 0 j15;
   j21 j22 0 0 0;
   j31 0 j31 0 0;
   0 0 0 0 j45;
   0 0 0 j54 0];

J2 = J^2;
% Jf = J*f0;

dj11=[-2*x0(1) 0 0 0 0];

dj11f=dj11*f0;

dJf=[dj11f 0 0 0 0;
         0 0 0 0 0;
         0 0 0 0 0;
         0 0 0 0 0;
         0 0 0 0 0];
     
dj11Jf=dj11*(J*f0);
dJJf=[dj11Jf 0 0 0 0;
         0 0 0 0 0;
         0 0 0 0 0;
         0 0 0 0 0;
         0 0 0 0 0];

%ddjf11 =    [-4*x0(1) + 8/3*x0(1)^3 + 2*x0(2) + 2*x0(3) - 2*x0(5), 2*x0(1), 2*x0(1), 0, -2*x0(1)];
ddff11=-2*f0(1)*f0(1);%ddjf11*f0;
ddJff=[ddff11 0 0 0 0;
            0 0 0 0 0;
            0 0 0 0 0;
            0 0 0 0 0;
            0 0 0 0 0];

ddf = (J2 + dJf );
%d3f = ddJ + 2*dJ*J + J*ddf;
%d3f = ddJff + 2*dJf*J + J*dJf + J^3;
%d3f = ddJff + dJf*J + 2*J*dJf + J^3;
%d3f = ddJff + dJJf + (2*dJf  + J2)*J + J*dJf;
d3f = ddJff + dJJf + 2*dJf*J + J*dJf + J^3;

Dnum4 = (I + h^2*(-1/4*J2 + 1/3*ddf))*h*f0;
Den4 = I - 0.5*h*J + h^2/6*ddf - h^3/24*d3f;


%Samelson inverse yeah:

% den = Den4*f0;
% 
% den_inv = den/(den'*den);
% 
% dx = (den_inv'*Dnum4)*f0;

%*********************

denCond = cond(Den4);

%[L,U] = lu(Den4);
%dx = U\(L\Dnum4); 

%matrix feature:

dx = Den4\Dnum4;

%***************

x = x0 + dx; %4 ord

%************ PREV ONE ***********
% Dnum3 = Dnum4 - h^4/24*d3f*f0;
% dx = U\(L\Dnum3);
%************ ALT ONE ************
Den3 = I - 0.5*h*J + h^2/6*ddf;
dx = Den3\Dnum4;
% 
x2 = x0 + dx; %3 ord

%************ BSRK ESTIMATOR *************
%x2 = BSRK23_step(fun,t0,h,x0);
%*********************************
J2=J;
err = x - x2;
%err = h^4/24*d3f*f0;

end