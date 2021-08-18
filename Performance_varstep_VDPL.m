function Performance_varstep_VDPL
%draw performane plot for VDPL test problem

 mu = 1000;
 function dx = F(t, x) 
     dx = x;
     dx(1) = x(2);
     dx(2) = mu*(1 - x(1)^2)*x(2) - x(1);
 end
 x0 = [2; 0]; %initial value
 tfin = 2500; %experiment time
 fun = @F; %poiner at F

h = 1e-5; %initial stepsize
t0 = 0;   %initial time

%time step limits
hmin = 1e-9;
hmax = 200;

Ntol = 7; %number of experimental points on "eps" axis
tolspan = logspace(-9,-2,Ntol);
Nexp = 7; %number of experiments for time measuring

%errors of solvers
% errRN = zeros(Nexp,Ntol);   %Rational numerical 4(3)
errLo = zeros(Nexp,Ntol);   %LobattoIIIC 4(3)
errRA = zeros(Nexp,Ntol);   %Rational analytical 4(3)
errTaylor = zeros(Nexp,Ntol); %Taylor 4(3)
errERK= zeros(Nexp,Ntol);   %ERK 4(3)
% errode15s = zeros(Nexp,Ntol); %ode15s

%latencies of solvers
% latRN = zeros(Nexp,Ntol);
latLo = zeros(Nexp,Ntol);
latRA = zeros(Nexp,Ntol);
latTaylor = zeros(Nexp,Ntol);
latERK= zeros(Nexp,Ntol);
% latode15s= zeros(Nexp,Ntol);

%start progress bar
hw = waitbar(0,'wait...');

%Get reference solution with DOPRI8(7)
abstol2 = 1e-14;
reltol2 = 1e-11;
hmin2 = 1e-11;

% [~, Xet] = DOPRI78(fun,t0,tfin,x0,h,hmin2,hmax,abstol2,reltol2);
% Xet = transpose(Xet);
options = odeset('AbsTol',abstol2,'RelTol',reltol2,'InitialStep',h,'MaxStep',hmax,'MaxOrder',5); 
[~, Xet] = ode15s(fun,[t0,tfin],x0,options);

for latctr = 1:Nexp
    for i = 1:Ntol
        reltol = tolspan(i);
        abstol = reltol*1e-5;
        
%         tic
%         %options = odeset('AbsTol',abstol,'RelTol',reltol,'InitialStep',h,'MaxStep',hmax,'MaxOrder',3);
%         options = odeset('AbsTol',abstol,'RelTol',reltol,'InitialStep',h,'MaxStep',hmax);
%         [~,Xode15s] = ode23s(fun,[t0,tfin],x0,options);
%         latode15s(latctr,i) = toc;
        
%         tic
%         [~,XP] = RN34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
%         latRN(latctr,i) = toc;
%         XP = transpose(XP); %to set dimension to the standard matlab style
        
        tic
        [~,XL] = Lobatto34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latLo(latctr,i) = toc;
        XL = transpose(XL);
        
        tic
        [~,XPa] = RA34_VDPL(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,mu);
        latRA(latctr,i) = toc;
        XPa = transpose(XPa); %to set dimension to the standard matlab style
        
        tic
        [~,XT] = Taylor34filtVDPL(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,mu);
        latTaylor(latctr,i) = toc;
        XT = transpose(XT); %to set dimension to the standard matlab style
        
        tic
        [~,XE] = ERK34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latERK(latctr,i) = toc;
        XE = transpose(XE); %to set dimension to the standard matlab style
        
        %find errors in the last point
%         errRN(latctr,i)    = -log(norm(Xet(end,:) - XP(end,:)));
        errLo(latctr,i)    = -log(norm(Xet(end,:) - XL(end,:)));
        errRA(latctr,i)    = -log(norm(Xet(end,:) - XPa(end,:)));
        errTaylor(latctr,i)= -log(norm(Xet(end,:) - XT(end,:)));
        errERK(latctr,i)   = -log(norm(Xet(end,:) - XE(end,:)));
%         errode15s(latctr,i)= -log(norm(Xet(end,:) - Xode15s(end,:)));
        
    end
    waitbar(latctr/Nexp,hw);
end
close(hw);

%averaging: errors
% errRN = median(errRN,1);
errLo = median(errLo,1);
errRA = median(errRA,1);
errTaylor = median(errTaylor,1);
errERK = median(errERK,1);
% errode15s = median(errode15s,1);

%averaging: latency
% latRN = median(latRN,1);
latLo = median(latLo,1);
latRA= median(latRA,1);
latTaylor = median(latTaylor,1);
latERK = median(latERK ,1);
% latode15s = median(latode15s ,1);

figure(1);
semilogy(errLo,latLo,'-s',errRA,latRA,'v-',errTaylor,latTaylor,'-o',errERK,latERK,'.-g');
ylabel('$T$','interpreter','latex');     
xlabel('$-log(err)$','interpreter','latex');
legend('LobattoIIIC4(3)','RA4(3)','Taylor4(3)','ERK4(3)','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

end

