function Test_varstep_HIRES
% Тест во временной области задачи Риккати
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
tfin = 100; %experiment time

h = 1e-3;
t0 = 0;
fun = @F;

%Varstep Solver
hmin = 1e-5;
hmax = 100;
abstol = 1e-10;
reltol = 1e-5;

disp('RA34');
tic
[tspanAnalyt,X1,~] = RA34_HIRES(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanAnalyt)


disp('LobattoIIIC 34');
tic
[tspanLo,XL] = Lobatto34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanLo)

disp('Taylor');
tic
[tspanTaylor,XTaylor] = Taylor34filtHIRES(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanTaylor)

disp('ERK34');
tic
[tspanERK,XERK] = ERK34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanERK)


abstol = 1e-16;
reltol = 1e-11;
hmin = 1e-11;

options = odeset('AbsTol',abstol,'RelTol',reltol);
disp('ode15s ref');
tic
[tet2, Xet2] = ode15s(fun,[t0,tfin],x0,options); Xet2 = Xet2';
toc


nvar = 1;
figure(3);

plot(tspanLo,XL(nvar,:),'--',tspanAnalyt,X1(nvar,:),'.-',tspanTaylor,XTaylor(nvar,:),'.-',tspanERK,XERK(nvar,:),'.-g',tet2,Xet2(nvar,:),'.-k');
legend('LobattoIIIC4(3)','RA4(3)','Taylor4(3)','ERK4(3)','DOPRI78 Reference','interpreter','latex');
ylabel('$y_1$','interpreter','latex');     
xlabel('$t$','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');


figure(2);
subplot(2,1,1);
plot(tspanLo,XL(1,:),'.-',tspanLo,XL(2,:),'.-',tspanLo,XL(3,:),'.-');
legend('$y_1$','$y_2$','$y_3$','interpreter','latex');
ylabel('$y$','interpreter','latex');     
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

