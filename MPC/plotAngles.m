%plot angles in time domain
function res = plotAngles(id, tout, xout, uout, ps)

    f = figure(id);
    f.Position = [600 0 1000 900];
    uout = [uout; uout(end,:)];

    %plot joint values
    subplot(3,1,1);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Angular position [rad]");
    title("Position of J1, J2");
    plot(tout,xout(:,1),"b-","LineWidth",1.5);
    plot(tout,xout(:,2),"r-","LineWidth",1.5);
    plot(ps.t,ps.xref(:,1), "b--");
    plot(ps.t,ps.xref(:,2), "r--");
    legend("$x_1$","$x_2$","Interpreter","latex");
    yticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]);
    yticklabels({'-\pi','-0.75\pi','-0.5\pi','-0.25\pi','0','0.25\pi','0.5\pi','0.75\pi','\pi'});

    %plot joint velocities
    subplot(3,1,2);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Angular velocity [rad/s]");
    title("Angular velocity of J1,J2");
    plot(tout,xout(:,3),"b-","LineWidth",1.5);
    plot(tout,xout(:,4),"r-","LineWidth",1.5);
    plot(ps.t,ps.xref(:,3), "b--");
    plot(ps.t,ps.xref(:,4), "r--");
    legend("$\dot{x_1}$", "$\dot{x_2}$","Interpreter","latex")

    %plot input
    subplot(3,1,3);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Input Value");
    title("Input");
    plot(tout,uout(:,1),"b-","LineWidth",1.5);
    plot(tout,uout(:,2),"r-","LineWidth",1.5);
    %ylim([-1.2 1.2])
    legend("$u_1$","$u_2$","Interpreter","latex")
    
    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);
    filename = sprintf("plots/state/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    saveas(id,filename);

    res = 0;
end

