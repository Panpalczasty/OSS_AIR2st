% Plot switch function with input
%
% Parameters:
%   id - id of figure
%   tout - time vector
%   uout - input to system (as array)
%   Huout - switch function (as array)
%   ps - model & system parameter struct
function res = plotSwitch(id, tout, uout, Huout, ps)
    
    uout = [uout; uout(end,:)];
    f = figure(id);
    f.Position = [0 0 600 400];

    subplot(2,1,1);
    plot(tout, uout(:,1),"b-", "LineWidth", 1.5);
    hold on;
    plot(tout, Huout(:, 1), "b--");

    grid;
    xlabel("Time [s]");
    ylabel("Value");
    title("Switch fcn - JT1");
    labels = ["$u_1$", "$\phi_1$"];
    legend(labels, "Interpreter","latex");

    subplot(2,1,2);
    plot(tout, uout(:,2), "r-", "LineWidth", 1.5);
    hold on;
    plot(tout, Huout(:, 2), "r--");

    grid;
    xlabel("Time [s]");
    ylabel("Value");
    title("Switch fcn - JT2");
    labels = ["$u_2$", "$\phi_2$"];
    legend(labels, "Interpreter","latex");

    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);
    filename = sprintf("plots/switch/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    saveas(id,filename);
 
    res = 0;
end
