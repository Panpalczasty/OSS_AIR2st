%gradient checker
function dif = checkGradient(ps, ep)
    
    nt = length(ps.u);

    gn = zeros(nt,1);

    for i = 1:nt
        dtau = zeros(nt,1);
        dtau(i) = ep;
        
        uright = ps.u + dtau;
        uleft = ps.u - dtau;
        
        [Qright, ~] = costFun(ps, uright);
        [Qleft, ~] = costFun(ps, uleft);
        
        gn(i) = (Qright - Qleft)/(2*ep);
       
    end

    [~, g] = costFun(ps, ps.taul);

    disp(["Ham","Num"]);
    disp([g gn]);

    dif = gn;
end