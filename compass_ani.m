% numerically simulating compass system with compass_Eu.m or compass_RK4.m
% first

% Animation of the chaotic compass system
% Ctrl+c to terminate, close the figure window after terminating to prevent
% figure error

if ismember({'Odiv','tot_ts','b2_peri_ts','tau','theta','b1','b2','wetau','delta'},who) == ones(1,9) % check the existence of needed variables
    t_ani_st = 100; % system time of starting animation
    figure ;
    nf = gcf ; % index of figure
    figure(nf)
    switch Odiv
        case {1}
            for m = (t_ani_st/tau):tot_ts
                if (mod(m,b2_peri_ts/50)==0)
                    compass(1.5*exp(1i*theta(m+2,1)),'b')
                    hold on
                    compass(b1/100.0,b2*cos(m*wetau+delta)/100.0,'k')
                    hold off
                    title(['simulation time:',num2str(m*tau,'%.2f'),'(T)'])
                    figure(nf)
                end
            end
        case {2}
            for m = (t_ani_st/tau):tot_ts
                if (mod(m,b2_peri_ts/50)==0)
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
            for m = (t_ani_st/tau):tot_ts
                if (mod(m,b2_peri_ts/50)==0)
                    compass(1.5*exp(1i*theta(m+2,1)),'b')
                    hold on
                    compass(b1/100.0,b2*cos(m*wetau+delta)/100.0,'k')
                    compass(1.5*exp(1i*theta(m+2,2)),'r')
                    compass(1.5*exp(1i*theta(m+2,3)),'g')
                    hold off
                    title(['simulation time:',num2str(m*tau,'%.2f'),'(T)'])
                    figure(nf)
                end
            end
    end
    
    hold off
end