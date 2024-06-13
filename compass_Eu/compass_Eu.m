% numerically simulating compass system with Euler method
% numerical method: Euler

clear all % clearing Workspace
tic % starting timer
tau = 0.0001 ; % length of time step
tot_time = 500 ; % total simulation time
tot_ts = tot_time/tau ; % total time steps
w_ext = 2*pi ; % angular velocity of driving field , T=2*pi/w_ext 
b2_peri_ts = round( (2*pi/w_ext)/tau) ; % times steps in a period of B field oscillation
gamma = 6.0 ; % dimensionless damping constant
b1 = 36.0 ; % dimensionless strength of static component of B field
b2 = 99.3 ; % dimensionless amplitude of oscillating component of B field. direction is prependicular to b1
Odiv = 2 ; % number of initial condition(s)
delta = 0 ; % initial phase of external B field in rad
%delta = rand(1,1)*2*pi ;
%delta_de = delta/pi*180 ;

% memory allocation
theta = zeros((tot_ts+2),Odiv) ;
w = zeros( (tot_ts+2),Odiv) ;
theta_n = zeros(floor(tot_time),Odiv) ;
theta_ave = zeros(floor(tot_time),Odiv) ;
w_n = zeros(floor(tot_time),Odiv) ;

% calculate some constant first
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;

for k = 1:Odiv
    
    % set initial condition
    w(1,k) = 0.0 + (k-1)*0.0  ;
    w(2,k) = 0.0 ;
    theta(1,k) = 0.0 + -1*(k-1)*(k-3)*0.25*2*pi ;
    theta(2,k) = theta(1,k) + w(1,k)*tau ;
    jj = 1 ; % indices of theta_n and <theta>_n
    for m = 1:tot_ts
        % numerical calculation with central difference Euler method 
        theta(m+2,k) = (theta(m,k)*(gata2-1.0) + theta(m+1,k)*2.0 + tausq*(-b1*sin(theta(m+1,k) )+b2*sin(pi/2-theta(m+1,k) )*cos(m*wetau+delta) ) )/(1+gata2) ;
        w(m+1,k) = (theta(m+2,k) - theta(m,k) )/(2*tau) ;
        if(mod(m,b2_peri_ts)==0)
            theta_n(jj,k) = theta(m+2,k) ;
            w_n(jj,k) = w(m+1,k) ;
            nc = floor((theta_n(jj,k) + pi)/2/pi) ; % translating (n*2pi) make angle value between (-pi)~(pi) 
            if( nc ~= 0)
                theta_n(jj,k) = theta_n(jj,k) - nc*2*pi ;
            end
            theta_ave(jj,k) = mean(theta( (m+2-b2_peri_ts+1):m+2,k) ) ;
            nc = floor((theta_ave(jj,k) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_ave(jj,k) = theta_ave(jj,k) - nc*2*pi ;
            end
            jj = jj+1 ;
        end
    end
    for m = 1:tot_ts
        nc = floor((theta(m+2,k) + pi)/2/pi) ;
        if( nc ~= 0)
            theta(m+2,k) = theta(m+2,k) - nc*2*pi ;
        end
    end

    w(tot_ts+2,k) = (theta(tot_ts+2,k)-theta(tot_ts+1,k) )/tau ;
    
   
    
    figure; plot(theta( ( (tot_time-400)/tau):100:(tot_time/tau),k)./(2*pi),w( (tot_time-400)/tau:100:(tot_time/tau),k)./(2*pi),'.','MarkerSize',8)
    xlabel('\theta')
    ylabel('\omega')
    title(['B_2=', sprintf('%.2f',b2),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
        
    %{
    % theta and w when time=integer
    hold on
    plot(theta_n(k,( (tot_time-800)):tot_time )./pi*180,w_n(k,(tot_time-800):tot_time )./pi*180,'r.','MarkerSize',2)
    hold off
    
    figure; plot(theta_n(k,( (tot_time-800)):tot_time )./pi*180,w_n(k,(tot_time-800):tot_time )./pi*180,'r.','MarkerSize',2)
    xlabel('\theta')
    ylabel('\omega')
    title(['B_2=', num2str(b2),', \theta_0=',num2str(theta(k,1)/pi*180 ),'\circ'])
    %}

    % Poincare map
    figure; plot(theta_n(tot_time-200:tot_time-1,k)./(2*pi),theta_n(tot_time-199:tot_time,k)./(2*pi),'.','MarkerSize',6)
    xlabel('\theta_n')
    ylabel('\theta_{n+1}')
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])

    figure; plot(theta_ave(tot_time-200:tot_time-1,k)./(2*pi),theta_ave(tot_time-199:tot_time,k)./(2*pi),'.','MarkerSize',6)
    xlabel('<\theta>_n')
    ylabel('<\theta>_{n+1}')
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])

    figure;plot( 0:tau:20,theta(tau/tau:20/tau+1,k)./(2*pi),'b')
    xlabel('time')
    ylabel('\theta')
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
    ylim([-0.5 0.5])

    figure;plot( (tot_time-20):tau:tot_time,theta( (tot_time-20)/tau:tot_time/tau,k)./(2*pi),'b')
    xlabel('time')
    ylabel('\theta')
    title(['B_2=',num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
    ylim([-0.5 0.5])

end

runtime = toc
