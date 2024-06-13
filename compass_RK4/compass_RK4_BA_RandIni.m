% Calculating Basin of Attraction(BA) with randomly choosing initial condition input by RK4 method

clear all
tic

tau = 0.0001 ;
tot_time_max = 300 ;
tot_ts_max = round(tot_time_max/tau) ; % ts:time step
w_ext = 2*pi ;
b2_per_ts = round(2*pi/w_ext/tau) ; % per:period
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 103.2 ;
n_ini = 2 ; % number of initial conditions
w_min = -1.75 ; % minimum of omega in (rad/(2*pi) )
w_max = 1.75 ; % maximum of omega in (rad/(2*pi) )
rng(33) % seed of Random Number Generator. artificial
theta0 = rand(n_ini,1)-0.5 ; % initial angle in (rad/(2*pi) )
w0 = rand(n_ini,1).*(w_max-w_min)+w_min ; % initial omega in (rad/(2*pi) )

n_baplot_a = 0 ; % count of 1st kind period attractor
n_baplot_b = 0 ; % count of 2nd kind period attractor
n_no_plot = 0 ; % count of not (1st or 2nd)
coss = 10 ; % condition of successive satisfying
cur_ssa = 0 ; % current successive satisfying for a
cur_ssb = 0 ; % current successive satisfying for b

% condition of identifying attractor. from other program. in rad
theta_ave_n3a = 0.0 ; % a condition. % here 3 stand for period-3
theta_ave_n5b = 0.0 ; % b condition. % here 5 stand for period-5

wetau = w_ext*tau ;

for k = 1:n_ini % loop change initial condition
    toc
    fprintf('Start %3.0f of %3.0f\n',k,n_ini)
    
    clear w theta theta_ave
    w = zeros( (tot_ts_max+1),1) ;
    theta = zeros( (tot_ts_max+1),1) ;
    theta_ave = zeros(tot_time_max,1) ;
    
    w(1) = w0(k)*2*pi ;
    theta(1) = theta0(k)*2*pi ;
    
    jj = 0 ; % count of theta_ave
    n = 0 ; % count of while loop
    leave_while_n = 0 ; % switch variable: 0:had not satisfied stable condition; 1:satisfied
    while (n <= tot_ts_max) && (leave_while_n == 0)
        n = n+1 ;
        rkdth1 = w(n) ;
        rkdw1 = -gamma*w(n) - b1*sin(theta(n) ) + b2*cos(theta(n) )*cos( n*wetau) ; 
        rkdth2 = rkdth1 + 0.5*tau*rkdw1 ;
        rkdw2 = -gamma*rkdth2 - b1*sin(theta(n) + 0.5*rkdth1*tau) + b2*cos(theta(n) + 0.5*rkdth1*tau)*cos( (n+0.5)*wetau) ; 
        rkdth3 = rkdth1 + 0.5*tau*rkdw2 ;
        rkdw3 = -gamma*rkdth3 - b1*sin(theta(n) + 0.5*rkdth2*tau) + b2*cos(theta(n) + 0.5*rkdth2*tau)*cos( (n+0.5)*wetau) ; 
        rkdth4 = rkdth1 + tau*rkdw3 ;
        rkdw4 = -gamma*rkdth4 - b1*sin(theta(n) + rkdth3*tau) + b2*cos(theta(n) + rkdth3*tau)*cos( (n+1)*wetau) ; 
        theta(n+1) = theta(n) + tau*(rkdth1 + 2*rkdth2 + 2*rkdth3 + rkdth4)/6.0 ;
        w(n+1) = w(n) + tau*(rkdw1 + 2*rkdw2 + 2*rkdw3 + rkdw4 )/6.0 ;
        if(mod(n,b2_per_ts)==0 )
            jj = jj+1 ;
            theta_ave(jj) = mean(theta(n+1-b2_per_ts:n+1) ) ;
            nc = floor((theta_ave(jj) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_ave(jj) = theta_ave(jj) - nc*2*pi ;
            end
            theta_ave(jj) = mean(theta((n+1-b2_per_ts+1):n+1) ) ;
            nc = floor((theta_ave(jj) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_ave(jj) = theta_ave(jj) - nc*2*pi ;
            end
            if (jj>20)
                if  (abs( mean(theta_ave(jj-2:jj) ) - theta_ave_n3a) < 10^-3) % condition for a. here use theta_ave=0. 10^-3 is artificial.
                    cur_ssa = cur_ssa + 1 ;
                    if (cur_ssa >= coss)
                        n_baplot_a = n_baplot_a + 1 ;
                        baplot_a(1,n_baplot_a) = theta(1) ;
                        baplot_a(2,n_baplot_a) = w(1) ;
                        leave_while_n = 1 ;
                    end
                elseif (abs( mean(theta_ave(jj-4:jj) ) - theta_ave_n5b) < 10^-3) % condition for b
                    cur_ssb = cur_ssb + 1 ;
                    if (cur_ssb >= coss)
                        n_baplot_b = n_baplot_b + 1 ;
                        baplot_b(1,n_baplot_b) = theta(1) ;
                        baplot_b(2,n_baplot_b) = w(1) ;
                        leave_while_n = 1 ;
                    end
                else
                    cur_ssa = 0 ;
                    cur_ssb = 0 ;
                end
            end
            
            
        end
    end
    
    % not (1st or 2nd)
    if(leave_while_n == 0)
        n_no_plot = n_no_plot + 1 ;
        baplot_no(1,n_no_plot) = theta(1) ;
        baplot_no(2,n_no_plot) = w(1) ;
    end

end
    

figure; hold on
if n_baplot_a~=0
    plot(baplot_a(1,:)/(2*pi),baplot_a(2,:)/(2*pi),'color','b','LineStyle','none','Marker','.','MarkerSize',6) ; % blue stands for 1st kind of periodic attractor
end
if n_baplot_b~=0
    plot(baplot_b(1,:)/(2*pi),baplot_b(2,:)/(2*pi),'color','r','LineStyle','none','Marker','.','MarkerSize',6) ; % red stands for 2nd kind of periodic attractor
end
if n_no_plot~=0
    plot(baplot_no(1,:)/(2*pi),baplot_no(2,:)/(2*pi),'color','k','LineStyle','none','Marker','.','MarkerSize',6) ; % black stands for not (1st or 2nd)
end
hold off
xlim([-0.5 0.5]);ylim([(w_min-0.1) (w_max+0.1)])
xlabel('\theta_0(X2\pi)')
ylabel('\omega_0(X2\pi/T)')
title(['B_2=', num2str(b2,'%.2f'), ' RK4'])

runtime = toc

