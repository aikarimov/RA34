function J = J_ord(fun, t0, x0, deltaX, ord)
%J_ORD finds partial derivatives with order
%ord = 2
%ord = 4

dim = length(x0);

if (ord == 2)
    dX = zeros(dim,1);
    dX(1) = deltaX;
    J = zeros(dim);
    for i = 1:dim
        J(:,i) = 0.5*(feval(fun,t0,x0 + dX) - feval(fun,t0,x0 - dX))/deltaX; %ord 2
        dX = circshift(dX,1);
    end
end

if (ord == 4)
    dX = zeros(dim,1);
    dX(1) = deltaX;
    J = zeros(dim);
    for i = 1:dim
        J(:,i) = (-1/12*feval(fun,t0,x0 + 2*dX) + 2/3*feval(fun,t0,x0 + dX)- 2/3*feval(fun,t0,x0 - dX) + 1/12*feval(fun,t0,x0 - 2*dX))/deltaX; %ord 4
        dX = circshift(dX,1);
    end
end

end