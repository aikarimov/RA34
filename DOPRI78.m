function [tspan, X] = DOPRI78(fun,t0,tfin,x0,h0,hmin,hmax,abstol, reltol)
% DOPRI78 Solve non-stiff differential equations with variable-step Dormand-Prince 8(7) method.
%     [TSPAN, X] = DOPRI78(FUN,T0,TFIN,X0,H0,HMIN,HMAX,ABSTOL, RELTOL)
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
    
    tol = max(reltol*norm(x),abstol);

    [xnew, errR] = DOPRI78_step(fun,t,h,x);
    errnorm = norm(errR);
    if (errnorm > tol) && (h > hmin)
        h = StepControl(h,hmin,hmax,tol,errnorm,7);
        rcount = rcount + 1;
    else
        t = t + h;
        h = StepControl(h,hmin,hmax,tol,errnorm,7);
        tspan(i) = t;
        x = xnew;
        X(:,i) = x;
        i = i + 1;
        
        if (errnorm > tol)
            warning('tolerance is not satisfied!');
        end
    end
end
h = tfin - t;
[x, ~] = DOPRI78_step(fun,t,h,x);
t = t + h;
tspan(i) = t;
X(:,i) = x;
%rcount;
end