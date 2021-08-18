function Performance_varstep_FNR
%draw performane plot for FNR test problem

global A w

a=[0.5,0.7,0.8,0.002,-1]';
A=0.083; w = 0.00003*2*3.14;

function dx = F(t, x)
    dx = x;
    dx(1)=x(1)-(x(1)*x(1)*x(1))/3-x(2)-x(3)+x(5);
    dx(2)=a(1)*(a(2)+x(1)-a(3)*x(2));
    dx(3)=a(4)*(a(5)-x(3)-x(1));
    dx(4)=A*x(5);
    dx(5)=-w^2*x(4);
end
x0 = [0;0;0;0;1];
tfin = 7e4;%1.4e4


%solver settings
h = 1e-2;
hmin = 1e-7;
hmax = 1e4;
t0 = 0;
fun = @F;


Ntol = 7; %points on "eps" axis
tolspan = logspace(-8,-4,Ntol);
Nexp = 7; %experiments

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

hw = waitbar(0,'wait...');

abstol2 = 1e-16;
reltol2 = 1e-12;
h0 = 1e-2;
tic
% [~, Xet] = DOPRI78(fun,t0,tfin,x0,h0,hmin,hmax,abstol2,reltol2);
% Xet = transpose(Xet);
options = odeset('AbsTol',abstol2,'RelTol',reltol2,'InitialStep',h,'MaxStep',hmax,'MaxOrder',5); 
[~, Xet] = ode15s(fun,[t0,tfin],x0,options);
toc

for latctr = 1:Nexp
    for i = 1:Ntol
        reltol = tolspan(i);
        abstol = reltol*1e-5;
        
        tic
        [~,XPa]  = RA34_FNR(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latRA(latctr,i) = toc
        XPa = transpose(XPa); %to set dimension to the standard matlab style
        
        tic
        [~,XT] = Taylor34filtFNR(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latTaylor(latctr,i) = toc
        XT = transpose(XT); %to set dimension to the standard matlab style
        
        tic
        [~,XE] = ERK34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latERK(latctr,i) = toc
        XE = transpose(XE); %to set dimension to the standard matlab style
        
        tic
        [~,XL] = Lobatto34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
        latLo(latctr,i) = toc
        XL = transpose(XL);
        
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

