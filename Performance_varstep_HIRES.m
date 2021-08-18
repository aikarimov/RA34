function Performance_varstep_HIRES
%draw performane plot for HIRES test problem

function dx = F(t,x)
    % F(t,x) is the derivative function for HIRES
    % t - time (not used)
    % x - state variable vector
   dx = x;
   dx(1) = -1.71*x(1) + 0.43*x(2) + 8.32*x(3) + 0.0007;
   dx(2) = 1.71*x(1) - 8.75*x(2);
   dx(3) = -10.03*x(3) + 0.43*x(4) + 0.035*x(5);
   dx(4) = 8.32*x(2) + 1.71*x(3) - 1.12*x(4);
   dx(5) = -1.745*x(5) + 0.43*x(6) + 0.43*x(7);
   dx(6) = -280*x(6)*x(8) + 0.69*x(4) + 1.71*x(5)  - 0.43*x(6) + 0.69*x(7);
   dx(7) = 280*x(6)*x(8) - 1.81*x(7);
   dx(8) = -280*x(6)*x(8) + 1.81*x(7);   
end
x0 = [1,0,0,0,0,0,0,0.0057]'; %initial value
tfin = 1000; %experiment time
fun = @F; %poiner at F

h = 1e-3; %initial stepsize
t0 = 0;   %initial time

%time step limits
hmin = 1e-9;
hmax = 200;

Ntol = 7;  %number of experimental points on "eps" axis
tolspan = logspace(-9,-2,Ntol);
Nexp = 7; %number of experiments for time measuring

%errors of solvers
errRN = zeros(Nexp,Ntol);   %Rational numerical 4(3)
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

%Get reference solution
abstol2 = 1e-15;
reltol2 = 1e-12;

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
        [~,XPa] = RA34_HIRES(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latRA(latctr,i) = toc;
        XPa = transpose(XPa); %to set dimension to the standard matlab style
        
        tic
        [~,XT] = Taylor34filtHIRES(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
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

