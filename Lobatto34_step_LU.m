function [x, err] = Lobatto34_step_LU(fun,t0,h,x0, tol)
%Lobatto IIIC RK method
global a c dim
a = [1/6, -1/3, 1/6;
    1/6, 5/12, -1/12;
    1/6, 2/3, 1/6];
c = [0, 0.5, 1];
b = [1/6, 2/3, 1/6];
b2 = [-1/2, 2, -1/2];

dim = length(x0);

dkvect = tol*100;

dx = tol*1;
kvect = zeros(3*dim,1);

Jval = J_ord(@(t,k)F(fun,t,k,h,x0),t0,kvect,dx,2); %numerical Jacobian

[L, U, P]  = lu(Jval);
%NEWTON ITERATIONS
while norm(dkvect) > tol %while change in k is larger then tolerance
    Fval = F(fun, t0, kvect, h, x0);
    y = L\(P*(-Fval));
    dkvect = U\y;
    kvect = kvect + dkvect;
end

x = x0;
x2 = x0;

for i=1:3
   kcur = kvect(1 + dim*(i-1):dim*i);
   x = x + h*b(i)*kcur;
   x2 = x2 + h*b2(i)*kcur;
end

err = x - x2;

end

function fnewt = F(fun, t0, kvect, h, x0)
% kvect = [ [k1] [k2] [k3]]';
% size(k) = [dim,3]
global a c dim

k = zeros(dim,3);
for i=1:3
   k(:,i) = kvect(1 + dim*(i-1):dim*i);
end

fval = zeros(3*dim,1);
ksum = zeros(dim,3);
for i = 1:3
    for j = 1:3
        ksum(:,i) = ksum(:,i) + a(i,j)*k(:,j);
    end
end

for i=1:3
   fval(1 + dim*(i-1):dim*i) = feval(fun,t0 + c(i)*h,x0 + h*ksum(:,i));
end

fnewt = kvect - fval;
end














