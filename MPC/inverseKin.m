%inverse kinematics - above or below
function [t, x] = inverseKin(pos, T, config, J1, J2)

    rad = sqrt(pos(:,1).^2 + pos(:,2).^2);
    b = (rad.^2 - J1.d^2 - J2.d^2)/(2*J1.d*J2.d);

    samples = length(pos(:,1));
    sfreq = samples(1)/T;
    t = linspace(0, T, samples);

    if config == "above"
       c = 1; 
    elseif config == "below"
       c = -1;
    end
 
    if any(abs(b)>1)  
        disp("Trajectory out of robot range!");
        x = zeros(samples,4);
    else
        alpha = atan2(-c*J2.d*sqrt(1-b.^2), J1.d + J2.d*b);
    
        x1 = alpha + atan2(pos(:,2), pos(:,1));
        x2 = c*acos(b); 
        x3 = [diff(x1)*sfreq; 0];
        x4 = [diff(x2)*sfreq; 0];    
    
        x3(end) = x3(end-1);
        x4(end) = x4(end-1);

        x = [x1 x2 x3 x4];
    end
end

%gradient based on fmincon - TODO
function [uout,Qstat] = fminconGrad(ps)

    Q=@(u) costFun(ps,u);
    
    options = optimoptions('fmincon');
    options.SpecifyObjectiveGradient = true;
    options.MaxFunctionEvaluations = 1e4;
    options.Display = 'iter';
    options.OptimalityTolerance = eps;
    options.StepTolerance = eps;
    options.Algorithm = 'interior-point';

    uopt = fmincon(Q,ps.u,ps.A,ps.b,[],[],[],[],[],options);
    
    Qstat = costFun(ps, uopt);

    %unravel
    uout(:,1) = uopt(1:end/2);
    uout(:,2) = uopt(end/2+1:end);
end