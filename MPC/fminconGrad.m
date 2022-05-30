%gradient based on fmincon
function [uout,Qstat] = fminconGrad(ps)

    Q=@(u) costFun(ps,u);
    
    options = optimoptions('fmincon');
    options.SpecifyObjectiveGradient = true;
    options.MaxFunctionEvaluations = 1e5;
    options.Display = 'off';
    options.OptimalityTolerance = eps;
    options.StepTolerance = eps;
    options.Algorithm = 'interior-point';

    uopt = fmincon(Q,ps.u,ps.A,ps.b,[],[],[],[],[],options);
    
    Qstat = costFun(ps, uopt);

    %unravel
    uout(:,1) = uopt(1:end/2);
    uout(:,2) = uopt(end/2+1:end);
end