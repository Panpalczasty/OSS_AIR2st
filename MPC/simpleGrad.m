%simple gradient descent algorithm
function [uopt, Qstat] = simpleGrad(ps, alpha, iter)
    
    figure(1)
    hold on;
    u1plot = plot(ps.posref(:,1), ps.posref(:,2), 'b--');
    u2plot = plot(ps.posref(:,1), ps.posref(:,2), 'k-',"LineWidth",1.5);
    legend("goal","trajectory");
    xlabel("time[s]");
    ylabel("input");
    grid;

    Qstat = zeros(iter,1);
    Qprev=0;
    Q = 10;
    i = 2;
    while (i <= iter) && (abs(Qprev-Q) >= 1e-7) 
        Qprev=Q;
        [Q, g] = costFun(ps, ps.u);
        ps.u = ps.u - g*alpha;
        
        %animate input
        [~, x] = qsolve45(ps.u, ps.T, ps);
        pos = simpleKin(x, ps.J1, ps.J2);
        set(u2plot,'XData',pos(:,1));
        set(u2plot,'YData',pos(:,2));
        drawnow
    
        Qstat(i) = Q;
        disp([i,Q])
        i = i+1;
    end
    uopt = ps.u;
end