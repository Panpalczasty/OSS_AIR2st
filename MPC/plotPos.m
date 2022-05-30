
%plot trajectory in 2D plane
function res = plotPos(id, pos, posref, ps)
    
    f = figure(id);
    f.Position = [0 600 600 400];


    plot(posref(:,1), posref(:,2),"b--");
    hold on;
    grid on;
    plot(pos(1,1),pos(1,2), "b*");
    plot(pos(:,1),pos(:,2), "k","LineWidth",1.2);
    xlabel("x[m]");
    ylabel("y[m]");
    xlim([-1.7 1.7])
    ylim([-1.7 1.7])
    title("Trajectory of robot");
    legend("Reference trajectory", "Start point", "Trajectory");

    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);
    filename = sprintf("plots/pos/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    %saveas(id,filename);

    res = 0;
end
