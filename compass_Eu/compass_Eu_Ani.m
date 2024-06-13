% Animation of the chaotic compass system, calculating and plotting at the
% same time

% Ctrl+c to terminate, close the figure window after terminating to prevent
% figure error

% numerical method: Euler
% black arrow indicates the direction and magnitude of magnetic field 
% other colors(1:blue, 2:red, 3:green) indicate the direction of magnetic dipole

clear all

tau = 0.0001 ; % length of time step
tot_time = 500 ; % total simulation time
tot_ts = round(tot_time/tau) ; % total time steps
w_ext = 2*pi ; % angular velocity of driving field , T=2*pi/w_ext 
b2_peri_ts = round( (2*pi/w_ext)/tau) ; % times steps in a period of the B field oscillation
gamma = 6.0 ; % dimensionless damping constant
b1 = 36.0 ; % dimensionless strength of fixed component of B field
b2 = 105.7 ; % dimensionless amplitude of oscillating component of B field. direction is prependicular to b1
Odiv = 1 ; % divide a circle to how many part. 1~3 is permissible in this program.
figure ;
nf = gcf ; % index of figure

delta = 0.0 ; % initial phase of external B field

% memory allocation
theta = zeros((tot_ts+2),Odiv) ;
w0 = zeros( 2,Odiv) ;
theta_n = zeros(floor(tot_time),Odiv) ;
theta_ave = zeros(floor(tot_time),Odiv) ;

% calculate some constant first
gata2 = gamma*tau/2.0 ;
tausq = tau*tau ;
wetau = w_ext*tau ;

% set initial conditions
if ( Odiv >= 1 )
    w0(1,1) = 0.0 ;
    w0(2,1) = 0.0 ;
    theta(1,1) = 0*2*pi ;
    theta(2,1) = theta(1,1) + w0(1,1)*tau ;
end
if (Odiv >= 2 )
    w0(1,2) = 0.0 ;
    w0(2,2) = 0.0 ;
    theta(1,2) = 0.278*2*pi ;
    theta(2,2) = theta(1,2) + w0(1,2)*tau ;
end
if (Odiv >= 3 )
    w0(1,3) = 0.0 ;
    w0(2,3) = 0.0 ;
    theta(1,3) = -0.3361*2*pi ;
    theta(2,3) = theta(1,3) + w0(1,3)*tau ;
end

% numerical calculation with central difference Euler method and plot
figure(nf)
switch Odiv
    case {1}
        for m = 1:tot_ts
            theta(m+2,:) = (theta(m,:)*(gata2-1.0) + theta(m+1,:)*2.0 + tausq*(-b1*sin(theta(m+1,:) )+b2*cos(theta(m+1,:) )*cos(m*wetau+delta) ) )/(1+gata2) ;

            if (mod(m,b2_peri_ts/50)==0 && m>0*b2_peri_ts)
                compass(1.5*exp(1i*theta(m+2,1)),'b')
                hold on
                compass(b1/100.0,b2*cos(m*wetau+delta)/100.0,'k')
                hold off
                title(['simulation time:',num2str(m*tau,'%.2f'),'(T)'])
                figure(nf)
                %pause(0.1)
            end

        end
    case {2}
        for m = 1:tot_ts
            theta(m+2,:) = (theta(m,:)*(gata2-1.0) + theta(m+1,:)*2.0 + tausq*(-b1*sin(theta(m+1,:) )+b2*cos(theta(m+1,:) )*cos(m*wetau+delta) ) )/(1+gata2) ;

            if (mod(m,b2_peri_ts/50)==0 && m>0*b2_peri_ts)
                compass(1.5*exp(1i*theta(m+2,1)),'b')
                hold on
                compass(b1/100.0,b2*cos(m*wetau+delta)/100.0,'k')
                compass(1.5*exp(1i*theta(m+2,2)),'r')
                hold off
                title(['simulation time:',num2str(m*tau,'%.2f'),'(T)'])
                figure(nf)
                %pause(0.1)
            end

        end
    case {3}
        for m = 1:tot_ts
            theta(m+2,:) = (theta(m,:)*(gata2-1.0) + theta(m+1,:)*2.0 + tausq*(-b1*sin(theta(m+1,:) )+b2*cos(theta(m+1,:) )*cos(m*wetau+delta) ) )/(1+gata2) ;

            if (mod(m,b2_peri_ts/50)==0 && m>0*b2_peri_ts)
                compass(1.5*exp(1i*theta(m+2,1)),'b')
                hold on
                compass(b1/100.0,b2*cos(m*wetau+delta)/100.0,'k')
                compass(1.5*exp(1i*theta(m+2,2)),'r')
                compass(1.5*exp(1i*theta(m+2,3)),'g')
                hold off
                title(['simulation time:',num2str(m*tau,'%.2f'),'(T)'])
                figure(nf)
                %pause(0.1)
            end
        end
end
hold off
