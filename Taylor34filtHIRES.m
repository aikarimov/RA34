function [tspan, X] = Taylor34filtVDPL(fun,t0,tfin,x0,h0,hmin,hmax,abstol, reltol, mu)
X = x0;
x = x0;
tspan = t0;
h = h0;
t = t0;

%deltaX = 0.001*abstol; %delta for Jacobi matrix eval

rcount = 0;

errnorm0 = [1 1];
%compute J - derivatives, Jacobian
%dim = length(x0);
%dX = zeros(dim,1);
%dX(1) = deltaX;
%J = J_ord(fun, t0, x0, deltaX, 2); %2 order Jacobian

i = 2; 
hp = h0;
% hp2 = hp;
% Jm1 = J; % J[n - 1]
% Jm2 = J; % J[n - 2]

NF = 3; %length of filtering sequence
Hspan = zeros(NF,1);
Espan = zeros(NF,1);

ord = 3;

while t < tfin - h
    %***** place Taylor methods here *****
    [xnew, errR] = Taylor34_step_1step_HIRES (fun,t,h,x);
    %*****
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
            h = StepControl(Hspan(1),hmin,hmax,tol,errnorm,ord);
        end
        
        rcount = rcount + 1;
        
        tspan(i) = t;
        x = xnew;
        X(:,i) = x;
        %denCond(i) = dCond;
        i = i + 1;
%         Jm2 = Jm1;
%         Jm1 = J;
%         J = J2;
        
        if (norm(errR) > tol)
            disp('warning: tolerance is not satisfied!');
        end
    end
    %shift back
    errnorm0(2) = errnorm(1);
    errnorm0(1) = errnorm;
   
end
hp2 = hp;
hp = h;
h = tfin - t;
[x, ~] =  Taylor34_step_1step_HIRES (fun,t,h,x);
t = t + h;
tspan(i) = t;
%denCond(i) = dCond;
X(:,i) = x;
%rcount
end