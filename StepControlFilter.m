function Hnew = StepControlFilter(Hspan,hmin,hmax,tol,Espan,k)
%filtering by IIR filters, based on E. Hairer's paper

    ksi = 0.9;
    gamma = 0.99;
    
    %Hnew = gamma*h0*(ksi*tol/err)^(1/((k+1)));
    h0 = Hspan(1);
    r = Hspan(1:end-1)./Hspan(2:end); %H[i]/H[i-1]
    eps = (ksi*tol./Espan); %(tol/e[i])^(1/(k+1))
    
    D = length(r);
    
    s = 2;
    
    if D == 2
        k = 16;
        b = [s*1/4/k, 1/2/k, (2 - s)*1/4/k]; 
        a = [1 0 0];
    end
    
    if D == 1
        b = [1/k, 0]; 
        a = [1 0]; 
    end
    
    rnew = (eps(1))^b(1);
    
    for i = 2:D+1
        rnew = rnew*(eps(i))^b(i)*(r(i-1))^(-a(i));
    end
    Hnew = gamma*h0*rnew;
        
    
    if Hnew < hmin
        Hnew = hmin;
    end
    if Hnew > hmax
        Hnew = hmax;
    end
end