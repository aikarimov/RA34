function [tspan, X, denCond] = RA34elwisefiltVDPL(fun,t0,tfin,x0,h0,hmin,hmax,abstol, reltol, mu)

X = x0;
x = x0;
tspan = t0;
h = h0;
t = t0;

rcount = 0;
denCond(1) = 1;
errnorm0 = [1 1];

i = 2; 
NF = 2; %length of filtering sequence
Hspan = zeros(NF,1);
Espan = zeros(NF,1);

ord = 3;

while t < tfin - h

    [xnew, errR, dCond] = RA34elwise_step_VDPL(fun,t,h,x,mu);

    
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
[x, ~, dCond] = RA34elwise_step_VDPL(fun,t,h,x,mu);
t = t + h;
tspan(i) = t;
denCond(i) = dCond;
X(:,i) = x;

end