function Test_varstep_FNR
% Variable-step FNR simulation
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
tfin = 7e4;

h = 1e-5;
t0 = 0;
fun = @F;

hmin = 1e-10;
hmax = 1000;

abstol = 1e-8;
reltol = 1e-5;

disp('RA34');
tic
[tspanAnalyt,X1,~] = RA34_FNR(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanAnalyt)

disp('LobattoIIIC 34');
tic
[tspanLo,XL] = Lobatto34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanLo)

disp('Taylor34');
tic
[tspanTaylor,XTaylor] = Taylor34filtFNR(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanTaylor)

disp('ERK34');
tic
[tspanERK,XERK] = ERK34(fun,t0,tfin,x0,h,hmin,hmax,abstol,reltol);
toc
length(tspanERK)

abstol = 1e-14;
reltol = 1e-11;

options = odeset('AbsTol',abstol,'RelTol',reltol);
disp('ode15s ref');
tic
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
plot(tet2,(Xet2(1,:)) ,'.-',tet2,(Xet2(2,:)) ,'.-',tet2,(Xet2(3,:)) ,'.-',tet2,(Xet2(4,:)/1000) ,'.-');
legend('$y_1$','$y_2$','$y_3$','$10^{-3}y_4$','interpreter','latex');

ylabel('$y$','interpreter','latex');     
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

