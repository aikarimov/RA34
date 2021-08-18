function [tspan, X] = Lobatto34(fun,t0,tfin,x0,h0,hmin,hmax,abstol, reltol)
% Lobatto34 solve differential equation with variable-step Explicit
% Runge-Kutta method, developed by Runge, Kutta, Dormand and Prince
%     [TSPAN, X] = Lobatto34(FUN,T0,TFIN,X0,H0,HMIN,HMAX,ABSTOL, RELTOL)
%   Input parameters:
%     FUN is a vector function FUN(t,x) = [f_1 ... f_M]'
%     T0 is starting time
%     TFIN is end time
%     X0 vector of initial conditions
%     HMIN is minimal stepsize
%     HMAX is maximal stepsize
%     ABSTOL is absolute error tolerance
%     RELTOL is relative error tolerance
%   Output parameters:
%     TSPAN - vector of time points, 1 x N
%     X is output array of solution, M x N
X = x0;
x = x0;
tspan = t0;
h = h0;
t = t0;
i = 2;

rcount = 0;

while t < tfin - h
    
    tol = max(reltol*norm(x), abstol);

    [xnew, errR] = Lobatto34_step_LU(fun,t,h,x, reltol);
    errnorm = norm(errR);
    if (errnorm > tol) && (h > hmin) %if too high error and we can make step lower, do not accept
        h = StepControl_2(h,hmin,hmax,tol,errnorm,2);
        rcount = rcount + 1;
    else
        t = t + h;
        h = StepControl_2(h,hmin,hmax,tol,errnorm,2);
        tspan(i) = t;
        x = xnew;
        X(:,i) = x;
        i = i + 1;
        
        if (errR > tol)
            disp('warning: tolerance is not satisfied!');
        end
    end
end
h = tfin - t;
[x, ~] = Lobatto34_step_LU(fun,t,h,x, reltol);
t = t + h;
tspan(i) = t;
X(:,i) = x;
%rcount
end