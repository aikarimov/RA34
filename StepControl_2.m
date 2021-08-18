function Hnew = StepControl_2(h0,hmin,hmax,tol,err,k)
    Hnew = 0.95*h0*(0.65*tol/err)^(1/(k+1));
    if Hnew < hmin
        Hnew = hmin;
    end
    if Hnew > hmax
        Hnew = hmax;
    end
end