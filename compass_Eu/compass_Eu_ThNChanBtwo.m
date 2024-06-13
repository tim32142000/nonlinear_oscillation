% plot theta(time=integer) versus different b2
% old technique: one b2 one plot
% numerical method: Euler

clear all
tau = 0.0001 ;
tot_time = 1500 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
tor_ext_ts = 2*pi/w_ext/tau ;
gamma = 6.0 ;
B1 = 36.0 ;
B2_min = 90.0 ;
B2_max = 120.0 ; % 117.77 one circle
B2_intv = 0.1 ;
B2_num = round( (B2_max-B2_min)/B2_intv + 1 ) ;
nop = 500 ; % number of points plotted per b2

gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;
B2 = B2_min ;
w(1) = 0.0 ;
w(2) = 0.0 ;
theta(1) = 85.0/180.0*pi ;
theta(2) = theta(1) + w(1)*tau ;

figure; 
hold on;
for n = 1:B2_num
    jj = 1 ;
    for m = 1:tot_ts
        theta(m+2) = (theta(m)*(gata2-1.0) + theta(m+1)*2.0 - tausq*(B1*sin(theta(m+1) )-B2*cos(theta(m+1) )*cos(m*wetau) ) )/(1+gata2) ;
        if(mod(m,tor_ext_ts)==0 )
            theta_n(jj) = theta(m+2) ;
            if( floor((theta_n(jj) + pi)/2/pi) ~= 0)
                nc = floor((theta_n(jj) + pi)/2/pi) ;
                theta_n(jj) = theta_n(jj) - nc*2*pi ;
            end
            jj = jj+1 ;
        end
    end
    plot(B2*ones(1,nop),theta_n(tot_time-nop+1:tot_time)./pi.*180,'LineStyle','none','Marker','.','MarkerEdgeColor','r','MarkerSize',6 )
    
    B2 = B2 + B2_intv ; % change b2
    
    % set next initial conditions of next parameter. here use final state of last parameter
    w(1) = w(length(w)-1) ;
    w(2) = w(length(w)) ;
    theta(1) = theta(length(theta)-1) ;
    theta(2) = theta(length(theta)) ;
    
end
hold off;
xlabel('B_2')
ylabel('\theta_n(\circ)')
