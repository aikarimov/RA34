function Test_varstep_VDPL
% Variable-step Van der Pol simulation
 mu = 1000;%150; %stiffness parameter
function dx = F(t, x) 
    dx = x;
    dx(1) = x(2);
    dx(2) = mu*(1 - x(1)^2)*x(2) - x(1);
end
x0 = [2; 0]; %initial value
tfin = 2000; %experiment time

%solver settings - for mu = 150
h = 1e-5;
t0 = 0;
fun = @F;

hmin = 1e-10;
hmax = 10;

abstol = 1e-11;
reltol = 1e-8;

disp('RA34_VDPL');
tic
[tspanAnalyt,X1,~] = RA34_VDPL(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,mu);
toc
length(tspanAnalyt)

disp('LobattoIIIC 34');
tic
[tspanLo,XL] = Lobatto34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanLo)

disp('Taylor34VDPL');
tic
[tspanTaylor,XTaylor] = Taylor34filtVDPL(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,mu);
%tspanTaylor = tspanAnalyt; XTaylor = X1;
toc
length(tspanTaylor)

disp('ERK34');
tic
[tspanERK,XERK] = ERK34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
%tspanERK = tspanAnalyt; XERK = X1;
toc
length(tspanERK)

abstol = 1e-14;
reltol = 1e-11;

options = odeset('AbsTol',abstol,'RelTol',reltol);
disp('ode15s ref');
tic
%[tet2, Xet2] = DOPRI78(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
[tet2, Xet2] = ode15s(fun,[t0,tfin],x0,options); Xet2 = Xet2';
toc

nvar = 1;
figure(3);

plot(tspanLo,XL(nvar,:),'--',...
    tspanAnalyt,X1(nvar,:),'.-',...
    tspanTaylor,XTaylor(nvar,:),'.-',...
    tspanERK, XERK(nvar,:),'.-g',...
    tet2,Xet2(nvar,:),'.-k');
legend('LobattoIIIC','RA34','Taylor','ERK34','ode15s Reference');


ylabel('$X$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
title('Time domain','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

figure(4);
subplot(2,1,1);
semilogy(tet2,abs(Xet2(1,:)) ,'.-',tet2,abs(Xet2(2,:)) ,'.-');
legend('$|y_1|$','$|y_2|$','interpreter','latex');

ylabel('$|y|$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
title('Solution','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

subplot(2,1,2);

semilogy(tspanLo(2:end),tspanLo(2:end) - tspanLo(1:end-1),'.-',...
         tspanAnalyt(2:end),tspanAnalyt(2:end) - tspanAnalyt(1:end-1),'.-',...
         tspanTaylor(2:end),tspanTaylor(2:end) - tspanTaylor(1:end-1),'.-',...
         tspanERK(2:end),tspanERK(2:end) - tspanERK(1:end-1),'.-g');
legend('LobattoIIIC4(3)','RA4(3)','Taylor4(3)','ERK4(3)','interpreter','latex');

width = 600;
height = 300;
set(gcf,'position',[100,100,width,height]);

title('Time step $h_k$','interpreter','latex');
ylabel('$h$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

disp('Lobatto error:');
disp(norm(Xet2(:,end) - XL(:,end)));

disp('Taylor error:');
disp(norm(Xet2(:,end) - XTaylor(:,end)));

disp('RA34 error:');
disp(norm(Xet2(:,end) - X1(:,end)));

disp('ERK34 error:');
disp(norm(Xet2(:,end) - XERK(:,end)));

end

