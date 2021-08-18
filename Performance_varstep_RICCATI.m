function Performance_varstep_RICCATI
%draw performane plot for RICCATI test problem
aparam = 10000;
function dx = F(t,X)
 %x0 = [0 0 1 0]';
    x1 = X(1); x2 = X(2); x3 = X(3); x4 = X(4);
    dx = zeros(4,1);
    dx(1) = aparam - x1^2 - x3*x2;
    dx(2) = - x2*(x1 + x4);
    dx(3) = - x3*(x1 + x4);
    dx(4) = aparam - x4^2 - x3*x2;
end
x0 = [0 0 1 0]'; %initial value
tfin = 10; %experiment time
fun = @F; %poiner at F

h = 1e-5; %initial stepsize
t0 = 0;   %initial time

%time step limits
hmin = 1e-9;
hmax = 200;

Ntol = 7; %number of experimental points on "eps" axis
tolspan = logspace(-8,-2,Ntol);
Nexp = 7; %number of experiments for time measuring


%errors of solvers
errLo = zeros(Nexp,Ntol);   %LobattoIIIC 4(3)
errRA = zeros(Nexp,Ntol);   %Rational analytical 4(3)
errTaylor = zeros(Nexp,Ntol); %Taylor 4(3)
errERK= zeros(Nexp,Ntol);   %ERK 4(3)

%latencies of solvers
latLo = zeros(Nexp,Ntol);
latRA = zeros(Nexp,Ntol);
latTaylor = zeros(Nexp,Ntol);
latERK= zeros(Nexp,Ntol);

%start progress bar
hw = waitbar(0,'wait...');

%Get reference solution with DOPRI8(7)
abstol2 = 1e-16;
reltol2 = 1e-12;
hmin2 = 1e-11;

% [~, Xet] = DOPRI78(fun,t0,tfin,x0,h,hmin2,hmax,abstol2,reltol2);
% Xet = transpose(Xet);

options = odeset('AbsTol',abstol2,'RelTol',reltol2,'InitialStep',h,'MaxStep',hmax,'MaxOrder',5); 
[~, Xet] = ode15s(fun,[t0,tfin],x0,options);

for latctr = 1:Nexp
    for i = 1:Ntol
        reltol = tolspan(i);
        abstol = reltol*1e-5;
              
        tic
        [~,XL] = Lobatto34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latLo(latctr,i) = toc;
        XL = transpose(XL);
        
        tic
        [~,XPa] = RA34_RICCATI(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,aparam);
        latRA(latctr,i) = toc;
        XPa = transpose(XPa); %to set dimension to the standard matlab style
        
        tic
        [~,XT] = Taylor34filtRICCATI(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,aparam);
        latTaylor(latctr,i) = toc;
        XT = transpose(XT); %to set dimension to the standard matlab style
        
        tic
        [~,XE] = ERK34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latERK(latctr,i) = toc;
        XE = transpose(XE); %to set dimension to the standard matlab style
        
        %find errors in the last point
        errLo(latctr,i)    = -log(norm(Xet(end,:) - XL(end,:)));
        errRA(latctr,i)    = -log(norm(Xet(end,:) - XPa(end,:)));
        errTaylor(latctr,i)= -log(norm(Xet(end,:) - XT(end,:)));
        errERK(latctr,i)   =-log(norm(Xet(end,:) - XE(end,:)));
        
    end
    waitbar(latctr/Nexp,hw);
end
close(hw);

%averaging: errors
errLo = median(errLo,1);
errRA = median(errRA,1);
errTaylor = median(errTaylor,1);
errERK = median(errERK,1);

%averaging: latency
latLo = median(latLo,1);
latRA= median(latRA,1);
latTaylor = median(latTaylor,1);
latERK = median(latERK ,1);

figure(1);
semilogy(errLo,latLo,'-s',errRA,latRA,'v-',errTaylor,latTaylor,'-o',errERK,latERK,'.-g');
ylabel('$T$','interpreter','latex');     
xlabel('$-log(err)$','interpreter','latex');
legend('LobattoIIIC4(3)','RA4(3)','Taylor4(3)','ERK4(3)','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

end

