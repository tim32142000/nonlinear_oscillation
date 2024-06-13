% plot theta(time=integer) versus different theta0
% old technique: one b2 one plot
% numerical method: Euler

clear all
tic
tau = 0.0001 ;
tot_time = 1000 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
b2_peri_ts = 2*pi/w_ext/tau ;
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 118.7 ;
Odiv = 90 ; % equally divide a circle to how many part
nop = 300 ; % number of points plotted per theta0

gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
for n = 1:1 % reserved for-loop for changing b2
    b2=b2 ;
    figure; hold on;
    for k = 1:Odiv
        fprintf('Start Odiv=%3.0f of %3.0f\n',k,Odiv)
        clear w theta theta_n
        w = zeros(1, (tot_ts+1)) ;
        theta = zeros(1, (tot_ts+1)) ;
        theta_n = zeros(1,tot_time) ;
        w(1) = 0.0 ;
        w(2) = 0.0 ;
        theta(1) = (k-1)/Odiv*2*pi ;
        theta(2) = theta(1) + w(1)*tau ;
        jj = 1 ;
        for m = 1:tot_ts
            theta(m+2) = (theta(m)*(gata2-1.0) + theta(m+1)*2.0 + tausq*(-b1*sin(theta(m+1) ) + b2*cos(theta(m+1) )*cos(m*wetau) ) )/(1+gata2) ;
            if(mod(m,b2_peri_ts)==0 )
                theta_n(jj) = theta(m+2) ;
                if( floor((theta_n(jj) + pi)/2/pi) ~= 0)
                    nc = floor((theta_n(jj) + pi)/2/pi) ;
                    theta_n(jj) = theta_n(jj) - nc*2*pi ;
                end
                jj = jj+1 ;
            end
        end
    plot(theta(1)/pi*180*ones(1,nop),theta_n(tot_time-nop+1:tot_time)./pi.*180,'LineStyle','none','Marker','.','MarkerSize',6 )
    %figure; plot(theta(k,(500/tau):(700/tau)-1),w(k,500/tau:(700/tau)-1),'.')
    end
    hold off;
    xlim([0 360]);ylim([-180 180])
    xlabel('\theta_0(\circ)')
    ylabel('\theta_n(\circ)')
    title(['B_2=', num2str(b2)])
end
runtime = toc