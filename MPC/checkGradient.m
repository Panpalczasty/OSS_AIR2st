%gradient checker
function dif = checkGradient(u, ps, ep)
    
    % first, ravel
    u_raveled = [u(:,1); u(:,2)];

    nu = length(u_raveled);

    gn = zeros(nu,1);

    for i = 1:nu
        du = zeros(nu,1);
        du(i) = ep;
        
        uright = u_raveled + du;
        uleft = u_raveled - du;
        
        [Qright, ~] = costFun(ps, uright);
        [Qleft, ~] = costFun(ps, uleft);
        
        gn(i) = (Qright - Qleft)/(2*ep);
       
    end

    [~, g] = costFun(ps, u_raveled);

    disp(["Ham","Num"]);
    disp([g gn]);

    dif = gn;

    figure(10)
    subplot(2,1,1);
    plot(g(1:end/2), 'r--');
    hold on;
    plot(gn(1:end/2), 'k', 'LineWidth', 1.2);
    grid on;
    legend("Hamiltonian", "Numeric");

    subplot(2,1,2);
    plot(g(end/2+1:end), 'r--');
    hold on;
    grid on;
    plot(gn(end/2+1:end), 'k', 'LineWidth', 1.2);
end