% Calculating Basin of Attraction(BA) with RK4 method

clear all
tic

tau = 0.0001 ;
tot_time_max = 1000 ;
tot_ts_max = round(tot_time_max/tau) ; % ts:time step
w_ext = 2*pi ;
b2_peri_ts = round(2*pi/w_ext/tau) ; % peri:period
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 95.0 ;
theta_min =  -0.5*2*pi; % in rad
theta_max =  0.5*2*pi ;
Odiv = 3 ; % equally divide a circle to how many part
w_min = -2*2*pi ; % in rad
w_max = 2*2*pi ;
wdiv = 3 ; % w equally deviding

theta0=linspace(theta_min,theta_max,Odiv) ;
w0=linspace(w_min,w_max,wdiv) ;
inicon=combvec(theta0,w0) ; % 2*N initial conditions matrix. 1st row:theta0, 2nd row:w0 

n_baplot_a = 0 ; % count of 1st kind period attractor
n_baplot_b = 0 ; % count of 2nd kind period attractor
n_no_plot = 0 ; % count of not (1st or 2nd) after tot_time_max
coss = 10 ; % condition of successive satisfying
cssa = 0 ; % current successive satisfying for a
cssb = 0 ; % current successive satisfying for b

% condition of identifying attractor. from other program. in rad
theta_ave1a = -0.4703420684169656 ; % here 1 stand for period-1
theta_ave1b = 0.4703420684169631 ; 

% calculating constant first
wetau = w_ext*tau ;

figure;
nf=gcf;

for nthw = 1:length(inicon(1,:) ) % loop change initial condition. nthw: index of initial conditions
    toc
    fprintf('Starting %.0f of %3.0f\n',nthw,length(inicon(1,:) ) )
    
    clear w theta theta_n theta_ave
    w = zeros(1, (tot_ts_max+1) ) ;
    theta = zeros(1, (tot_ts_max+1) ) ;
    theta_n = zeros(1,tot_time_max) ;
    theta_ave = zeros(1,tot_time_max) ;
    
    theta(1) = inicon(1,nthw) ;
    w(1) = inicon(2,nthw) ;
    
    jj = 0 ; % count of theta_n
    m = 0 ; % count of while loop
    leave_while_m = 0 ; % switch variable: 0:had not satisfied all condition; 1:satisfied
    while (m < tot_ts_max) && (leave_while_m == 0)
        m = m+1 ;
        rkdth1 = w(m) ;
        rkdw1 = -gamma*w(m) - b1*sin(theta(m) ) + b2*cos(theta(m) )*cos( m*wetau) ; 
        rkdth2 = rkdth1 + 0.5*tau*rkdw1 ;
        rkdw2 = -gamma*rkdth2 - b1*sin(theta(m) + 0.5*rkdth1*tau) + b2*cos(theta(m) + 0.5*rkdth1*tau)*cos( (m+0.5)*wetau) ; 
        rkdth3 = rkdth1 + 0.5*tau*rkdw2 ;
        rkdw3 = -gamma*rkdth3 - b1*sin(theta(m) + 0.5*rkdth2*tau) + b2*cos(theta(m) + 0.5*rkdth2*tau)*cos( (m+0.5)*wetau) ; 
        rkdth4 = rkdth1 + tau*rkdw3 ;
        rkdw4 = -gamma*rkdth4 - b1*sin(theta(m) + rkdth3*tau) + b2*cos(theta(m) + rkdth3*tau)*cos( (m+1)*wetau) ; 
        theta(m+1) = theta(m) + tau*(rkdth1 + 2*rkdth2 + 2*rkdth3 + rkdth4)/6.0 ;
        w(m+1) = w(m) + tau*(rkdw1 + 2*rkdw2 + 2*rkdw3 + rkdw4 )/6.0 ;
        if(mod(m,b2_peri_ts)==0 )
            jj = jj+1 ;
            theta_n(jj) = theta(m+1) ;
            nc = floor((theta_n(jj) + pi)/(2*pi) ) ;
            if( nc ~= 0)
                theta_n(jj) = theta_n(jj) - nc*2*pi ;
            end
            theta_ave(jj) = mean(theta((m+1-b2_peri_ts+1):m+1) ) ;
            nc = floor((theta_ave(jj) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_ave(jj) = theta_ave(jj) - nc*2*pi ;
            end
            if (jj>20) % after some evolving time, starting identifying
                if  (abs(theta_ave(jj) - theta_ave1a) < 10^-8) % condition for period1a. here use theta_ave. 10^-8 is artificial.
                    cssa = cssa + 1 ;
                    if (cssa >= coss)
                        n_baplot_a = n_baplot_a + 1 ;
                        baplot_a(1,n_baplot_a) = theta(1) ;
                        baplot_a(2,n_baplot_a) = w(1) ;
                        leave_while_m = 1 ;
                    end
                elseif (abs(theta_ave(jj) - theta_ave1b) < 10^-8) % 10^-8 is artificial
                    cssb = cssb + 1 ;
                    if (cssb >= coss)
                        n_baplot_b = n_baplot_b + 1 ;
                        baplot_b(1,n_baplot_b) = theta(1) ;
                        baplot_b(2,n_baplot_b) = w(1) ;
                        leave_while_m = 1 ;
                    end
                else
                    cssa = 0 ;
                    cssb = 0 ;
                end
            end
        end
    end
    
    % not (1st or 2nd)
    if(leave_while_m == 0)
        n_no_plot = n_no_plot + 1 ;
        baplot_no(1,n_no_plot) = theta(1) ;
        baplot_no(2,n_no_plot) = w(1) ;
    end
    

    if n_baplot_a > 0
        plot(baplot_a(1,:)/(2*pi),baplot_a(2,:)/(2*pi),'color','b','LineStyle','none','Marker','.','MarkerSize',6) ; % blue stands for 1st kind of periodic attractor
        hold on
    end
    if n_baplot_b > 0
        plot(baplot_b(1,:)/(2*pi),baplot_b(2,:)/(2*pi),'color','r','LineStyle','none','Marker','.','MarkerSize',6) ; % red stands for 2nd kind of periodic attractor
        hold on
    end
    if n_no_plot > 0
        plot(baplot_no(1,:)/(2*pi),baplot_no(2,:)/(2*pi),'color','k','LineStyle','none','Marker','.','MarkerSize',6) ; % black stands for not (1st or 2nd)
        hold on
    end
    if (Odiv>1)&&(wdiv>1)
        xlim([(theta_min-(theta_max-theta_min)/(Odiv-1) )/(2*pi) (theta_max+(theta_max-theta_min)/(Odiv-1) )/(2*pi)]);ylim([w_min/(2*pi)-(w_max-w_min)/( (wdiv-1)*2*pi) w_max/(2*pi)+(w_max-w_min)/( (wdiv-1)*2*pi)])
    end
    xlabel('\theta_0(\times2\pi)')
    ylabel('\omega_0(\times2\pi/T)')
    title(['B_2=', sprintf('%.2f',b2), ' RK4'])
    hold off
    figure(nf)
end

hold off

runtime = toc
