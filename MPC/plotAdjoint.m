%plot time travel parameters
function res = plotAdjoint(id, tout, pout, ps, fname)
    
    f = figure(id);
    f.Position = [0 0 500 500];

    %plot position adjoints
    subplot(2,1,1);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Weight component");
    title("Adjoint of position of J1, J2");
    plot(tout,pout(:,1),"b--","LineWidth",1.5);
    plot(tout,pout(:,2),"r--","LineWidth",1.5);
    legend("$\psi_1$","$\psi_2$","Interpreter","latex")

    %plot joint velocities
    subplot(2,1,2);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Weight component");
    title("Adjoint of velocity of J1,J2");
    plot(tout,pout(:,3),"b--","LineWidth",1.5);
    plot(tout,pout(:,4),"r--","LineWidth",1.5);
    legend("$\psi_3$", "$\psi_4$", "Interpreter","latex")

    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);

    if nargin > 4
        filename = sprintf("plots/adj/%s.png",fname);
    else
        filename = sprintf("plots/adj/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    end
    
    saveas(id,filename);

    res = 0;
end
