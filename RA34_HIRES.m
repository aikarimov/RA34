function [tspan, X, denCond] = RA34_HIRES(fun,t0,tfin,x0,h0,hmin,hmax,abstol, reltol)
%     [TSPAN, X, denCond] = RA34_HIRES(FUN,T0,TFIN,X0,H0,HMIN,HMAX,ABSTOL, RELTOL)
%   Input parameters:
%     FUN is a vector function FUN(t,x) = [f_1 ... f_M]' IMPORTANT: while
%       it accepts two parameters, it is assumed that the problem is
%       autonomus! FUN = f(x(t)), not f(t,x(t))!
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
%     denCond is denominator condition number

X = x0;
x = x0;
tspan = t0;
h = h0;
t = t0;


rcount = 0;

%compute J - derivatives, Jacobian
denCond(1) = 1;
%UNCOMMENT  FOR CHECKING NUMERICAL VS. ANALYTICAL
%deltaX = 0.001*abstol; %delta for Jacobi matrix eval
%J = J_ord(fun, t0, x0, deltaX, 2); %2 order Jacobian

i = 2; 

NF = 3; %length of filtering sequence
Hspan = zeros(NF,1);
Espan = zeros(NF,1);

ord = 3;

while t < tfin - h
    [xnew, errR, ~, dCond] = RA34_step_HIRES(fun,t,h,x);
    
    tol = max(reltol*norm(xnew),abstol);
    
    errnorm = norm(errR);
    if (errnorm > tol) && (h > hmin) %if too high error and we can make step lower, do not accept
        %shift step and error
        Espan(2:end) = Espan(1:end-1);
        Espan(1) = errnorm;
        
        Hspan(2:end) = Hspan(1:end-1);
        Hspan(1) = h;
        
        if rcount >= NF
            h = StepControlFilter(Hspan,hmin,hmax,tol,Espan,ord);
        else
            h = StepControl(h,hmin,hmax,tol,errnorm,ord);
        end
        rcount = rcount + 1;
    else
        %make a step
        t = t + h;
        
        %shift error and step
        Espan(2:end) = Espan(1:end-1);
        Espan(1) = errnorm;
        Hspan(2:end) = Hspan(1:end-1);
        Hspan(1) = h;
       
        if rcount >= NF
            h = StepControlFilter(Hspan,hmin,hmax,tol,Espan,ord);
        else
            h = StepControl(h0,hmin,hmax,tol,errnorm,ord);
        end
        
        rcount = rcount + 1;
        
        tspan(i) = t;
        x = xnew;
        X(:,i) = x;
        denCond(i) = dCond;
        i = i + 1;
        if (norm(errR) > tol)
            warning('tolerance is not satisfied!');
        end
    end
end
h = tfin - t;
[x, ~, ~,dCond] = RA34_step_HIRES (fun,t,h,x);
t = t + h;
tspan(i) = t;
denCond(i) = dCond;
X(:,i) = x;
end