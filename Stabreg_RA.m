% STABILITY REGIONS OF RATIONAL APPROXIMATION METHODS

%domain in a complex plane
sp = -15:0.1:5;
wp = -10:0.1:10;
[ss, ww] = meshgrid(sp,wp);

N = length(sp);
Rh = zeros(N,N,6);
%doulbe loop through the complex plane
for i = 1:N
    for j = 1:N
        s = ss(i,j);
        w = ww(i,j);
        
        z = s + 1i*w; % z - is the stability function parameter
        
        %RA methods stability regions, order 2-7
        Rh(i,j,1) = abs(1 + (z)/(1 - 0.5*z));
        Rh(i,j,2) = abs(1 + (z - 0.25*z^3 + 1/3*z^3)/(1 - 0.5*z + 1/6*z^2 ));
        Rh(i,j,3) = abs(1 + (z - 0.25*z^3 + 1/3*z^3)/(1 - 0.5*z + 1/6*z^2 - 1/24*z^3));
        Rh(i,j,4) = abs((1 + 0.5*z + 1/6*z^2 + 1/24*z^3 + 1/120*z^4 + 1/360*z^5)/(1 - 0.5*z + 1/6*z^2 - 1/24*z^3 + 1/120*z^4));
        Rh(i,j,5) = abs((1 + 0.5*z + 1/6*z^2 + 1/24*z^3 + 1/120*z^4 + 1/720*z^5)/(1 - 0.5*z + 1/6*z^2 - 1/24*z^3 + 1/120*z^4 - 1/720*z^5));
        Rh(i,j,6) = abs((1 + 0.5*z + 1/6*z^2 + 1/24*z^3 + 1/120*z^4 + 1/720*z^5 + 1/5040*z^6 + 2/40320*z^7)/(1 - 0.5*z + 1/6*z^2 - 1/24*z^3 + 1/120*z^4 - 1/720*z^5 + 1/5040*z^6));
    end
end

figure(1);
for i = 1:6 %loop through stability regions
    
    subplot(2,3,i);
    [c, hc] = contourf(ss,ww,1 - Rh(:,:,i),[1 0:0.01:1]); % filled shapes: lightest = 0, darkest = 1
    colormap('gray');
    set(hc,'edgecolor','none');

    %set description elements on the plot
    xlabel('$Re$','interpreter','latex');
    ylabel('$Im$','interpreter','latex');
    title(['$p = ',num2str(i+1),'$'],'interpreter','latex');
    set(gca,'ticklabelinterpreter','latex');
    axis equal
end