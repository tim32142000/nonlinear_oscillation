% plot 'time averages of theta in one driving period' versus different b2
% numerical method: Euler


clear all
tic
tau = 0.0001 ;
tot_time = 1500 ;
tot_ts = tot_time/tau ;
w_ext = 2*pi ;
tor_ext_ts = 2*pi/w_ext/tau ;
gamma = 6.0 ;
b1 = 36.0 ;
b2_min = 88.5 ;
b2_max = 88.7 ; % 117.77 over one circle
b2_intv = 0.1 ;
b2_num = round( (b2_max-b2_min)/b2_intv + 1 ) ;
nop = 300 ; % number of points plotted per b2

gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;

% set initial condition and memory allocation
b2 = b2_min ;
w(1) = 0.0 ;
w(2) = 0.0 ;
theta = zeros(tot_ts+2,1) ;
theta(1) = 0.0/180.0*pi ;
theta(2) = theta(1) + w(1)*tau ;
theta_ave = zeros(1,tot_time) ;
thavepl = zeros(2,b2_num*nop) ; % theta average plot
theta0 = theta(1) ;
w0 = w(1) ;

for n = 1:b2_num
    fprintf('Starting b2=%6.2f, %5.0f of %5.0f\n',b2,n,b2_num)
    toc
    jj = 0 ;
    thavepl(1,(n-1)*nop+1:n*nop) = b2 ;
    for m = 1:tot_ts
        theta(m+2) = (theta(m)*(gata2-1.0) + theta(m+1)*2.0 - tausq*(b1*sin(theta(m+1) )-b2*cos(theta(m+1) )*cos(m*wetau) ) )/(1+gata2) ;
        if(mod(m,tor_ext_ts)==0 )
            jj = jj+1 ;
            theta_ave(jj) = mean(theta(m+2-tor_ext_ts+1:m+2) ) ;
            nc = floor( (theta_ave(jj) + pi)/(2*pi) ) ;
            if( nc ~= 0)
                theta_ave(jj) = theta_ave(jj) - nc*2*pi ;
            end
            if(jj>(tot_time-nop))
                thavepl(2,nop*(n-1)+(jj-(tot_time-nop) ) ) = theta_ave(jj) ;
            end
        end
    end

    b2 = b2 + b2_intv ; % change b2
    
    % set next initial conditions of next parameter. here use final state of last parameter
    w(1) = w(length(w)-1) ;
    w(2) = w(length(w)) ;
    theta(1) = theta(length(theta)-1) ;
    theta(2) = theta(length(theta)) ;
    
end
figure; plot(thavepl(1,:),thavepl(2,:)./pi*180,'LineStyle','none','Marker','.','MarkerEdgeColor','r','MarkerSize',2)
xlim([b2_min-b2_intv b2_max+b2_intv])
title(['\theta_0=',num2str(theta0),', \omega_0=',num2str(w0),', # of T=',num2str(tot_time),', Euler'])
toc
