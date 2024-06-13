% caculating and plotting theta_n versus delta(initial phase of (co)sinusoidally oscillatiing magnetic field)
% old technique: one delta one plot
% numerical method: RK4

clear all
tic
tau = 0.0001 ;
tot_time = 800 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
tor_ext_ts = 2*pi/w_ext/tau ;
gamma = 6.0 ;
B1 = 36.0 ;
B2 = 95.00 ;
Odiv = 2 ; % equally divide a circle to how many part
nop = 300 ; % number of points plotted per delta

wetau = w_ext*tau ;
for n = 1:1 % reserved for-loop for changing b2
    B2=B2 ;
    figure; hold on;
    for k = 1:Odiv
        clear w theta theta_n
        w = zeros(1, (tot_ts+1)) ;
        theta = zeros(1, (tot_ts+1)) ;
        theta_n = zeros(1,tot_time) ;
        w(1) = 0.0 ;
        theta(1) = 143.0/180.0*pi ;
        delta = -pi + (k-1)/Odiv*2*pi ;
        jj = 1 ;
        for m = 1:tot_ts
            rkdth1 = w(m) ;
            rkdw1 = -gamma*w(m) - B1*sin(theta(m) ) - B2*cos(theta(m) )*cos( m*wetau+delta) ;
            rkdth2 = rkdth1 + 0.5*tau*rkdw1 ;
            rkdw2 = -gamma*rkdth2 - B1*sin(theta(m) + 0.5*rkdth1*tau) - B2*cos(theta(m) + 0.5*rkdth1*tau)*cos( (m+0.5)*wetau+delta) ;
            rkdth3 = rkdth1 + 0.5*tau*rkdw2 ;
            rkdw3 = -gamma*rkdth3 - B1*sin(theta(m) + 0.5*rkdth2*tau) - B2*cos(theta(m) + 0.5*rkdth2*tau)*cos( (m+0.5)*wetau+delta) ;
            rkdth4 = rkdth1 + tau*rkdw3 ;
            rkdw4 = -gamma*rkdth4 - B1*sin(theta(m) + rkdth3*tau) - B2*cos(theta(m) + rkdth3*tau)*cos( (m+1)*wetau+delta) ;
            theta(m+1) = theta(m) + tau*(rkdth1 + 2*rkdth2 + 2*rkdth3 + rkdth4)/6.0 ;
            w(m+1) = w(m) + tau*(rkdw1 + 2*rkdw2 + 2*rkdw3 + rkdw4 )/6.0 ;
            if(mod(m,tor_ext_ts)==0 )
                theta_n(jj) = theta(m+1) ;
                if( floor((theta_n(jj) + pi)/2/pi) ~= 0)
                    nc = floor((theta_n(jj) + pi)/2/pi) ;
                    theta_n(jj) = theta_n(jj) - nc*2*pi ;
                end
                jj = jj+1 ;
            end
        end
    plot(delta/pi*180*ones(1,nop),theta_n(tot_time-nop+1:tot_time)./pi.*180,'LineStyle','none','Marker','.','MarkerSize',6 )
    end
    hold off;
    xlim([-180-360/Odiv 180]);ylim([-180 180])
    xlabel('\delta(\circ)')
    ylabel('\theta_n(\circ)')
    title(['B_2=', num2str(B2)])
end
toc