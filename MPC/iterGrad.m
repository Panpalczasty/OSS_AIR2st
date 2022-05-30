function [uopt, Qstat] = iterGrad(ps, alpha, iter)
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
        %animate input
        %go forward in time and observe model
        [t, x] = qsolve45(ps.u, ps.T, ps);
        %get time machine
        pf = [- ps.Wstat*(x(end,1:end-1) - ps.xf)'; 0; 0];
        %go back in time and de-observe model
        p = qantisolve45(pf, t, x, ps.u, ps);
        Hu = getHuvec(t, x, p, ps.u, ps);
        ps.u = Hu + ps.R*ps.u;
        
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