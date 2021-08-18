function Test_varstep_RICCATI
% Тест во временной области задачи Риккати
global aparam;
aparam = 10000;
% ************ Riccati ************
function dx = F(t,X)
 %x0 = [0 0 1 0]';
    x1 = X(1); x2 = X(2); x3 = X(3); x4 = X(4);
    dx = zeros(4,1);
    dx(1) =aparam - x1^2 - x3*x2;
    dx(2) = - x2*(x1 + x4);
    dx(3) = - x3*(x1 + x4);
    dx(4) =aparam - x4^2 - x3*x2;
end
x0 = [0 0 1 0]';
tfin = 10;

h = 1e-3;
t0 = 0;
fun = @F;

%Varstep Solver
hmin = 1e-5;
hmax = 10;
abstol = 1e-10;
reltol = 1e-5;

disp('RA34');
tic
[tspanAnalyt,X1,~] = RA34_RICCATI(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,aparam);
toc
length(tspanAnalyt)


disp('LobattoIIIC 34');
tic
[tspanLo,XL] = Lobatto34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanLo)

disp('Taylor');
tic
[tspanTaylor,XTaylor] = Taylor34filtRICCATI(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol,aparam);
toc
length(tspanTaylor)

disp('ERK34');
tic
[tspanERK,XERK] = ERK34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanERK)


abstol = 1e-16;
reltol = 1e-11;

options = odeset('AbsTol',abstol,'RelTol',reltol);
disp('ode15s ref');
tic
[tet2, Xet2] = ode15s(fun,[t0,tfin],x0,options); Xet2 = Xet2';
toc


nvar = 1;
figure(3);

plot(tspanLo,XL(nvar,:),'--',tspanAnalyt,X1(nvar,:),'.-',tspanTaylor,XTaylor(nvar,:),'.-',tspanERK,XERK(nvar,:),'.-g',tet2,Xet2(nvar,:),'.-k');
legend('LobattoIIIC4(3)','RA4(3)','Taylor4(3)','ERK4(3)','ode15s Reference','interpreter','latex');
ylabel('$y_1$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');


figure(2);
subplot(2,1,1);
plot(tspanLo,abs(XL(1,:)),'.-',tspanLo,abs(XL(3,:)),'.-');
legend('$|y_1|$','$|y_3|$','interpreter','latex');
ylabel('$|y|$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
title('Solution','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

subplot(2,1,2);
semilogy(tspanLo(2:end),tspanLo(2:end) - tspanLo(1:end-1),'.-',...
    tspanAnalyt(2:end),tspanAnalyt(2:end) - tspanAnalyt(1:end-1),'.-',...
    tspanTaylor(2:end),tspanTaylor(2:end) - tspanTaylor(1:end-1),'.-',...
    tspanERK(2:end),tspanERK(2:end) - tspanERK(1:end-1),'.-g');
legend('LobattoIIIC4(3)','RA4(3)','Taylor4(3)','ERK34','interpreter','latex');

title('Time step $h_k$','interpreter','latex');
ylabel('$h$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

width = 600;
height = 300;
set(gcf,'position',[100,100,width,height]);

disp('Lobatto error:');
disp(norm(Xet2(:,end) - XL(:,end)));
disp('RA error:');
disp(norm(Xet2(:,end) - X1(:,end)));
disp('Taylor error:');
disp(norm(Xet2(:,end) - XTaylor(:,end)));
disp('ERK34 error:');
disp(norm(Xet2(:,end) - XERK(:,end)));
end

