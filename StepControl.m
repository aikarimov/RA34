function Hnew = StepControl(h0,hmin,hmax,tol,err,k)
    Hnew = 0.85*h0*(0.5*tol/err)^(1/((k+1)));
    if Hnew < hmin
        Hnew = hmin;
    end
    if Hnew > hmax
        Hnew = hmax;
    end
end