function fig = plotError(i, t, x, xref, ps, fname)
    
    err = abs(x(:,1:4)-xref);
    fig = figure(i);
    
    % plot angular position error
    subplot(3,1,1);
    plot(t, err(:,1), 'b-', 'LineWidth', 1.5);
    hold on;
    grid;
    plot(t, err(:,2), 'r-', 'LineWidth', 1.5);
    
    xlabel("Time [s]");
    ylabel("Angle Error [rad]");
    title("State error - angle");
    labels = ["$$|x_1 - x_{ref1}|$$", "$$|x_2 - x_{ref2}|$$"];
    legend(labels, Interpreter="latex");
    
    
    % plot angular velocity error
    subplot(3,1,2);
    plot(t, err(:,3), 'b-', 'LineWidth', 1.5);
    hold on;
    grid;
    plot(t, err(:,4), 'r-', 'LineWidth', 1.5);
    
    xlabel("Time [s]");
    ylabel("Angular Velocity Error [rad/s]");
    title("State error - angular velocity");
    labels = ["$$|\dot{x}_1 - \dot{x}_{ref1}|$$", "$$|\dot{x}_2 - \dot{x}_{ref2}|$$"];
    legend(labels, Interpreter="latex");
    
    % get cartesian pos
    posref = simpleKin(xref, ps.J1, ps.J2);
    pos = simpleKin(x, ps.J1, ps.J2);
    poserr = sqrt((pos(:,1)-posref(:,1)).^2 + (pos(:,2)-posref(:,2)).^2);
    
    % plot cartesian error
    subplot(3,1,3)
    
    plot(t, poserr, 'k-', 'LineWidth', 1.5);
    hold on;
    grid;
    
    xlabel("Time [s]");
    ylabel("Error [m]");
    title("Trajectory error");

    if nargin > 5
        filename = sprintf("plots/err/%s.png",fname);
    else
        filename = sprintf("plots/err/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    end

end

