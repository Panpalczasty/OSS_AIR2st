%adjoint variables checker
function dif = checkAdjoint(p, ps, step)
    
    xleft = ps.x0*ones(1,4) - eye(4,4)*step;
    xright = ps.x0*ones(1,4) + eye(4,4)*step;
    
    Qleft = zeros(4,1);
    Qright = zeros(4,1);

    for i = 1:length(ps.x0)
        ps.x0 = xleft(:,i);
        Qleft(i) = costFun(ps, ps.taul);

        ps.x0 = xright(:,i);
        Qright(i) = costFun(ps, ps.taul);
    end
    
    dQ = (Qright - Qleft)/(2*step);
    dif = p(1,:)' + dQ;
    disp(["Right", "Left", "Deriv", "Adj(0)", "Diff"])
    disp([Qright Qleft dQ -p(1,:)' dif])
end