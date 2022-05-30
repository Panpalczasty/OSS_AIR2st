%concatenate tau arrays - ?
function [tau,u] = tauconcat(taul, u0)
    
    T = taul(end);
    taul = taul(1:end-1);
    tau1 = taul(1:end/2);
    tau2 = taul(end/2+1:end);

    ntau = size(tau1,1) + size(tau2,1) + 1;
    nu = length(u0);
    
    tau = zeros(ntau, 1);
    u = zeros(ntau+1, nu);
    u(1,:) = u0;

    for i = 2:ntau
        [tau1cand, i1] = min(tau1);
        [tau2cand, i2]= min(tau2);
        tau(i-1) = min(tau1cand, tau2cand);

        if(tau1cand < tau2cand)
            tau1(i1) = inf;
            u(i, :) = [-u(i-1,1) u(i-1,2)];
        else
            tau2(i2) = inf;
            u(i, :) = [u(i-1,1) -u(i-1,2)];
        end
    end

    tau(end) = T;
    u(end,:) = u(end-1,:);
end