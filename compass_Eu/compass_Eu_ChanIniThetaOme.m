% Calculating theta-w Basin of Attraction(BA) with Euler method

clear all
tic
tau = 0.0001 ;
tot_time_max = 500 ; % max simulation time for per initial condition
tot_ts_max = tot_time_max/tau ; % max time step
w_ext = 2*pi ;
b2_peri_ts = round(2*pi/w_ext/tau) ; % time steps in a period of B
gamma = 6.0 ;
b1 = 36.0 ;
b2 = 99.3 ;
Odiv = 10 ; % equally divide a circle to how many part

% range of initial conditions in (degree/T)
w_min_de = -650 ;
w_max_de = 650 ;
w_intv_de = 1300 ; % intv:interval

% (degree/T) -> (rad/T)
w_min = w_min_de/180.0*pi ;
w_max = w_max_de/180.0*pi ;
w_intv = w_intv_de/180.0*pi ;

n_baplot_a = 0 ; % count of 1st kind period attractor
n_baplot_b = 0 ; % count of 2nd kind period attractor
coss = 10 ; % condition of successive satisfying
cssa = 0 ; % current successive satisfying for a
cssb = 0 ; % current successive satisfying for b

% condition of identifying attractor. from other program. in rad
theta_n2a = [-0.505826174094475 -0.050935426265787] ; % here 2 stand for period-2
theta_n2b = [-1.395779873278010 -1.228883903717313] ; 


% calculating some constant first
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;

    
%figure; hold on;
for p = w_min_de:w_intv_de:w_max_de % loop change w_0
    toc
    fprintf('Start w=%3.0f\n',p) % print processing w_0 value
    for k = 1:Odiv % loop change theta_0
        fprintf('Start Odiv=%3.0f of %3.0f\n',k,Odiv) % print processing theta_0 value
        
        % memory allocation
        clear w theta theta_n theta_ave
        w = zeros(1, 2) ;
        theta = zeros(1, (tot_ts_max+1)) ;
        theta_n = zeros(1,tot_time_max) ;
        theta_ave = zeros(1,tot_time_max) ;
        
        % set initial condition
        w(1) = p/180.0*pi ;
        theta(1) = -pi + (k-1)/Odiv*2*pi ;
        theta(2) = theta(1) + w(1)*tau ;            
           
        jj = 0 ; % count of theta_n
        m = 0 ; % count of while loop
        leave_while_m = 0 ; % switch variable: 0:has not satisfied stable condition; 1:has satisfied
        while m <= tot_ts_max && leave_while_m == 0
            m = m+1 ;
            theta(m+2) = (theta(m)*(gata2-1.0) + theta(m+1)*2.0 + tausq*(-b1*sin(theta(m+1) )+b2*cos(theta(m+1) )*cos(m*wetau) ) )/(1+gata2) ; % Euler solve ODE
            if(mod(m,b2_peri_ts)==0 ) % time=integer
                jj = jj+1 ;
                theta_n(jj) = theta(m+2) ;
                
                % translating (n*2pi) make angle value between (-pi)~(pi) 
                nc = floor((theta_n(jj) + pi)/2/pi) ;
                if( nc ~= 0)
                    theta_n(jj) = theta_n(jj) - nc*2*pi ;
                end
                if jj>= 2
                    if ( ( (abs(theta_n(jj) - theta_n2a(2) ) < 10^-4) && (abs(theta_n(jj-1) - theta_n2a(1) ) < 10^-4) )||( (abs(theta_n(jj-1) - theta_n2a(2) ) < 10^-4) && (abs(theta_n(jj) - theta_n2a(1) ) < 10^-4) ) ) % two conditions for period-2, and here use theta_n. 10^-4 is artificial.
                        cssa = cssa + 1 ;
                        if (cssa >= coss)
                            n_baplot_a = n_baplot_a + 1 ;
                            baplot_a(1,n_baplot_a) = theta(1) ;
                            baplot_a(2,n_baplot_a) = w(1) ;
                            leave_while_m = 1 ;
                        end
                    elseif ( ( (abs(theta_n(jj) - theta_n2b(2) ) < 10^-4) && (abs(theta_n(jj-1) - theta_n2b(1) ) < 10^-4) )||( (abs(theta_n(jj-1) - theta_n2b(2) ) < 10^-4) && (abs(theta_n(jj) - theta_n2b(1) ) < 10^-4) ) ) 
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
            
    end
end
    
figure; hold on
if n_baplot_a > 0
    plot(baplot_a(1,:)/pi*180,baplot_a(2,:)/pi*180,'color','b','LineStyle','none','Marker','.','MarkerSize',6) ;
end
if n_baplot_b >0
    plot(baplot_b(1,:)/pi*180,baplot_b(2,:)/pi*180,'color','r','LineStyle','none','Marker','.','MarkerSize',6) ;
end
hold off
xlim([-180-360/Odiv 180]);ylim([(w_min_de-w_intv_de) (w_max_de+w_intv_de)])
xlabel('\theta_0(\circ)')
ylabel('\omega_0(\circ/time)')
title(['B_2=', sprintf('%.2f',b2),' Euler'])
        


runtime = toc
