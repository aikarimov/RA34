function Test_varstep_VDPL_elwise

mu = 1000;
    function dx = F(t, x)
        dx = x;
        dx(1) = x(2);
        dx(2) = mu*(1 - x(1)^2)*x(2) - x(1);
    end
x0 = [2; 0]; %initial value
tfin = 2000; %experiment time

h = 1e-3;
t0 = 0;
fun = @F;

%Varstep Solver

hmin = 1e-10;
hmax = 10;

abstol = 1e-11;
reltol = 1e-8;

disp('Pade34');
tic
[tspanELW,X,denCond] = RA34elwisefiltVDPL(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,mu);

toc
length(tspanELW)

disp('RA34_VDPL');
tic
[tspanAnalyt,X1,denCond1] = RA34_VDPL(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,mu);
toc
length(tspanAnalyt)


disp('Taylor34VDPL');
tic
[tspanTaylor,XTaylor] = Taylor34filtVDPL(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,mu);
toc
length(tspanTaylor)



abstol = 1e-13;
reltol = 1e-10;
hmin = 1e-11;

disp('dopri78 ref');
tic
[tet2, Xet2] = DOPRI78(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc

nvar = 1;

figure(4);
subplot(2,1,1);
plot(tspanELW,X(nvar,:) ,'.-',tspanAnalyt,X1(nvar,:) ,'.-',tet2,Xet2(nvar,:),'.-');
legend('RA4(3) element-wise','RA4(3)','Dormand-Prince 8(7)','interpreter','latex');

ylabel('$y_1$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
title('Solution','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

subplot(2,1,2);

semilogy(tspanELW(2:end),tspanELW(2:end) - tspanELW(1:end-1),'.-',...
    tspanAnalyt(2:end),tspanAnalyt(2:end) - tspanAnalyt(1:end-1),'.-',...
    tspanTaylor(2:end),tspanTaylor(2:end) - tspanTaylor(1:end-1),'.-');
legend('RA4(3) element-wise','RA4(3)','Taylor4(3)','interpreter','latex');
title('Time step $h_k$','interpreter','latex');
ylabel('$h$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');


width = 600;
height = 300;
set(gcf,'position',[100,100,width,height]);

end

