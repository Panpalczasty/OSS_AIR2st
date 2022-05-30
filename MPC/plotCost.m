%plot cost
function res = plotCost(id, Qstat, ps)
    
    f = figure(id);
    f.Position = [0 0 600 400];

    semilogy(Qstat,"k-","LineWidth",2);
    grid;
    xlabel("Iteration");
    ylabel("Cost value");
    title("Cost function stats");

    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);
    filename = sprintf("plots/cost/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    saveas(id,filename);
 
    res = 0;
end