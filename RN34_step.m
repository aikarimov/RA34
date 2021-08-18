function [x, err, J2, denCond] = RN34_step(fun,t0,h,x0,deltaX,J)
% RN 4(3) step mehtod implementation for autonomous problems
% first step
% fun - function
% t0 - time
% h = h[0] - current step
% deltaX - scalar delta for Jacobian comutation
% J   = J[ 0] Jacobian matrix

dim = length(x0);
I = eye(dim);

f0 = feval(fun,t0,x0);

H = 1e-1*h;
%H = 1e-8;
%Compute by ERK4 x[n+1] ord4 and x[n+1/2] ord4

[xp1, ~] = ERK34_step(fun,t0,H,x0);
[xp12, ~] = ERK34_step(fun,t0, 0.5*H,x0);

% H = h;
% [xp1, ~] = DOPRI78_step(fun,t0,H,x0);
% [xp12, ~] = DOPRI78_step(fun,t0, 0.5*H,x0);


%Compute by explicit RK ord 2

% [xm12, ~] = ERK34_step(fun,t0,-0.5*H,x0);
% [xp12, ~] = ERK34_step(fun,t0, 0.5*H,x0);
% 
% [xm1, ~] = ERK34_step(fun,t0,-H,x0);
% [xp1, ~] = ERK34_step(fun,t0, H,x0);

% 
% k1 = f0;
% x12 = x0 - 0.25*H*k1;
% k2 = feval(fun,t0 - 0.25*H,x12);
% xm12 = x0 - 0.5*H*k2; %2 ord
% 
% k1 = f0;
% x14 = x0 + 0.25*H*k1;
% k2 = feval(fun,t0 + 0.25*H,x14);
% xp12 = x0 + 0.5*H*k2; %2 ord
% 
% Jm12  = J_ord(fun,t0 - 0.5*H,xm12,deltaX,2);
% Jp12  = J_ord(fun,t0 + 0.5*H,xp12,deltaX,2);
% dJ = 0.25/H*(-Jm12 + Jp12);%coefs 1/2, 1/2*h
% ddJ = 0.25/H^2*(Jm12 - 2*J + Jp12);


% 
% k1 = f0;
% x12 = x0 + 0.5*H*k1;
% k2 = feval(fun,t0 + 0.5*H,x12);
% xp1 = x0 + H*k2; %2 ord
% 
% k1 = f0;
% x14 = x0 + 0.25*H*k1;
% k2 = feval(fun,t0 + 0.25*H,x14);
% xp12 = x0 + 0.5*H*k2; %2 ord

% 
Jp1   = J_ord(fun,t0 +     H,xp1,deltaX,2);
Jp12  = J_ord(fun,t0 + 0.5*H,xp12,deltaX,2);
% Jm1   = J_ord(fun,t0     - H,xm1,deltaX,4);
% Jm12  = J_ord(fun,t0 - 0.5*H,xm12,deltaX,4);

%ord 2
dJ = 0.5/H*(-3/2*J + 2*Jp12 - 0.5*Jp1);
ddJ = 0.25/H^2*(J - 2*Jp12 + Jp1);

%ord4
% dJ =   0.5/H  *( 1/12*Jm1 - 2/3*Jm12         + 2/3*Jp12 - 1/12*Jp1);
% ddJ = 0.25/H^2*(-1/12*Jm1 + 4/3*Jm12 - 5/2*J + 4/3*Jp12 - 1/12*Jp1);


df = J;
ddf = (J^2 + dJ);
%d3f = ddJ + 2*dJ*J + J*ddf;
d3f = ddJ + dJ*J + 2*J*dJ + J^3;%********** NEW VER ***********

Dnum4 = (I + h^2*(-1/4*df^2 + 1/3*ddf))*h*f0;
Den4 = I - 0.5*h*df + h^2/6*ddf - h^3/24*d3f;

denCond = cond(Den4);

%[L,U] = lu(Den4);
%dx = U\(L\Dnum4);
dx = Den4\Dnum4;
x = x0 + dx; %4 ord

%************ PREV ONE ***********

%Dnum3 = Dnum4 - h^4*d3f*f0/24;
%dx = U\(L\Dnum3);

%************ ALT ONE ************
 Den3 = I - 0.5*h*J + h^2/6*ddf;% + h^3/6*d3f;
 dx = Den3\Dnum4;
% 
%*********************************

x2 = x0 + dx; %3 ord

%************ BSRK ESTIMATOR *************
%x2 = BSRK23_step(fun,t0,h,x0);
%*********************************

J2  = J_ord(fun,t0 + h,x,deltaX,4);

err = x - x2;
%err = h^4/24*d3f*f0;
%err = max(abs(h^3*ddf*f0/6) , abs(h^4*d3f*f0/24)); 
end