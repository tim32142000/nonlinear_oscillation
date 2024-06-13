% solve compass system's ODEs with numerical method: RK4

clear all
tic
tau = 0.0001 ; % time step
tot_time = 1000 ; % system evolve time
tot_ts = round(tot_time/tau) ; % total time step
w_ext = 2*pi ; % driving angular frequency
b2_peri_ts = round( (2*pi/w_ext)/tau) ; % number of time steps in a period of external oscillatory magnetic field
gamma = 6.0 ; % damping
b1 = 36.0 ; % strength of fixed part
b2 = 106.6; % amplitude of oscillatory part
Odiv = 2 ; % number of initial condition(s)
delta = 0 ; % initial phase of external B field in rad

alpha = pi/2 ; % angle from B_1 direction to B_2 ccw-direction

theta = zeros( (tot_ts+1),Odiv) ;
w = zeros( (tot_ts+1),Odiv) ;
theta_n = zeros(floor(tot_time),Odiv) ;
theta_ave = zeros(floor(tot_time),Odiv) ;
w_n = zeros(floor(tot_time),Odiv) ;

wetau = w_ext*tau ;
for k = 1:Odiv
    theta(1,k) = (k-2)*(k-3)/2*0*2*pi + -1*(k-1)*(k-3)*0.25*2*pi + (k-1)/2.0*(k-2)*-0.3361*2*pi ;
    w(1,k) = (k-2)*(k-3)/2*0*2*pi + -1*(k-1)*(k-3)*0*2*pi ;
    jj = 0 ; % indices of theta_n and <theta>_n
    for m = 1:tot_ts % use rk4 to solve ODEs
        rkdth1 = w(m,k) ;
        rkdw1 = -gamma*rkdth1 - b1*sin(theta(m,k) ) + b2*sin(alpha-theta(m,k) )*cos( m*wetau+delta) ; 
        rkdth2 = rkdth1 + 0.5*tau*rkdw1 ;
        rkdw2 = -gamma*rkdth2 - b1*sin(theta(m,k) + 0.5*rkdth1*tau) + b2*sin(alpha-(theta(m,k) + 0.5*rkdth1*tau) )*cos( (m+0.5)*wetau+delta) ;
        rkdth3 = rkdth1 + 0.5*tau*rkdw2 ;
        rkdw3 = -gamma*rkdth3 - b1*sin(theta(m,k) + 0.5*rkdth2*tau) + b2*sin(alpha-(theta(m,k) + 0.5*rkdth2*tau) )*cos( (m+0.5)*wetau+delta) ; 
        rkdth4 = rkdth1 + tau*rkdw3 ;
        rkdw4 = -gamma*rkdth4 - b1*sin(theta(m,k) + rkdth3*tau) + b2*sin(alpha-(theta(m,k) + rkdth3*tau) )*cos( (m+1)*wetau+delta) ; 
        theta(m+1,k) = theta(m,k) + tau*(rkdth1 + 2*rkdth2 + 2*rkdth3 + rkdth4)/6.0 ;
        w(m+1,k) = w(m,k) + tau*(rkdw1 + 2*rkdw2 + 2*rkdw3 + rkdw4 )/6.0 ;
        if(mod(m,b2_peri_ts)==0) % t=integer
            jj = jj + 1 ;
            theta_n(jj,k) = theta(m+1,k) ; % save Poincare section theta_n
            w_n(jj,k) = w(m+1,k) ;
            nc = floor((theta_n(jj,k) + pi)/2/pi) ;
            if( nc ~= 0) % translate angle by(2n*pi) to [-pi,pi]
                theta_n(jj,k) = theta_n(jj,k) - nc*2*pi ;
            end
            theta_ave(jj,k) = mean(theta( (m+1-b2_peri_ts+1):m+1,k) ) ; % save averaged theta over one period
            nc = floor((theta_ave(jj,k) + pi)/2/pi) ;
            if( nc ~= 0)
                theta_ave(jj,k) = theta_ave(jj,k) - nc*2*pi ;
            end
        end
    end
    for m = 1:tot_ts 
        nc = floor( (theta(m+1,k) + pi)/2/pi) ;
        if( nc ~= 0)
                theta(m+1,k) = theta(m+1,k) - nc*2*pi ;
        end
    end

    % theta-w phase portrait
    figure; plot(theta( ( (tot_time-100)/tau):100:(tot_time/tau),k)./(2*pi),w( (tot_time-100)/tau:100:(tot_time/tau),k)./(2*pi),'.','MarkerSize',2)
    %{
    hold on
    plot(theta_n(tot_time-300:tot_time,k)./(2*pi),w_n(tot_time-300:tot_time,k)./(2*pi),'ro','MarkerSize',6)
    hold off
    %}
    xlabel('$\theta$(rad/$2\pi$)','interpreter','latex','fontsize',30)
    ylabel('$\dot \theta$(rad/$2\pi/T_d)$','interpreter','latex','fontsize',30)
    %title(['$b_2=', sprintf('%.2f',b2),', \theta_0=',num2str(theta(1,k)/(2*pi) ),'$'],'interpreter','latex','fontsize',16)
   
    %{
    % theta_{n+1} versus theta_n
    figure; plot(theta_n(tot_time-300:tot_time-1,k)./(2*pi),theta_n(tot_time-299:tot_time,k)./(2*pi),'.','MarkerSize',6)
    xlabel('\theta_n')
    ylabel('\theta_{n+1}')
    title(['B_2=', num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])

	% <theta>_{n+1} versus <theta>_n
    figure; plot(theta_ave(tot_time-300:tot_time-1,k)./(2*pi),theta_ave(tot_time-299:tot_time,k)./(2*pi),'.','MarkerSize',6)
    xlabel('<\theta>_n')
    ylabel('<\theta>_{n+1}')
    title(['B_2=',num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])

	% theta-t earlier part
    figure;plot( 0:tau:20,theta(tau/tau:20/tau+1,k)./(2*pi),'b')
    xlabel('time(T)')
    ylabel('\theta')
    title(['B_2=',num2str(b2,'%.2f'),', \theta_0=',num2str(theta(1,k)/(2*pi) )])
    ylim([-0.5 0.5])
    %}
    
	% theta-t with (t mod 1)
    %{
    figure;plot(mod( (tot_time-20):tau:tot_time,1),theta( (tot_time-20)/tau:tot_time/tau,k)./(2*pi),'b.','MarkerSize',2)
    xlabel('time(T)')
    ylabel('\theta')
    title(['B_2=', sprintf('%.2f',b2),', \theta_0=',num2str(theta(1,k)/pi*180.0 )])
    ylim([-0.5 0.5])
    hold on
    plot( (0:tau:1),atan(b2*cos(wetau.*(0:10000) )./b1)./(2*pi),'k')
    hold off
    %}
    
	% theta-t
    figure;plot( (tot_time-20):100*tau:tot_time,theta( (tot_time-20)/tau:100:tot_time/tau,k)./(2*pi),'b.','MarkerSize',2)
    xlabel('Time($T_d$)','interpreter','latex','fontsize',30)
    ylabel('$\theta$(rad/$2 \pi$)','interpreter','latex','fontsize',30)
    %title(['$b_2=', sprintf('%.2f',b2),', \theta_0=',num2str(theta(1,k)/(2*pi) ),'$'],'interpreter','latex','fontsize',16)
    ylim([-0.5 0.5])
    
end

runtime = toc
