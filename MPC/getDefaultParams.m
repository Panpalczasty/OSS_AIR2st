%declare modelfunction [t, x] = inverseKin(pos, T, config, J1, J2)
function params = getDefaultParams(params)

    %define params for joints
    
    %Define joint 1
    J1.I =  1.00; %inertia [kgm^2]
    J1.m =  10.0; %mass [kg]
    J1.d =  1.00; %length [m]
    J1.xc=  0.50; %x center coord [m]
    J1.yc=  0.00; %y center coord [m]
    J1.s =  0.10; %max static friction torque [Nm]
    J1.f =  0.01; %dynamic friction coeff [W-1]
    J1.um=  -1.00; %min steering torque [Nm]
    J1.uM=  1.00; %max steering torque [Nm]
    
    %Define joint 2
    J2.I =  1.00; %inertia [kgm^2]
    J2.m =  5.0; %mass [kg]
    J2.d =  0.70; %length [m]
    J2.xc=  0.35; %x center coord [m]
    J2.yc=  0.00; %y center coord [m]
    J2.s =  0.10; %max static friction torque [Nm]
    J2.f =  0.01; %dynamic friction coeff [W-1]
    J2.um=  -1.00; %min steering torque [Nm]
    J2.uM=  1.00; %max steering torque [Nm]
    
    %structify params
    params.J1 = J1;
    params.J2 = J2;
 
    %default trajectory
    t = linspace(0, params.T, params.T*params.sfreq);
    pos(:,1) = 1.7*ones(params.T*params.sfreq, 1);
    pos(:,2) = zeros(params.T*params.sfreq, 1);
    [~, params.xref] = inverseKin(pos, ...
                                  params.T, ...
                                  'above', ...
                                  params.J1, ...
                                  params.J2 ...
    )

    %initial and final state
    params.x0 = params.xref(1, :);
    params.xf = params.xref(end, :);

    %initial input - const
    params.u = zeros(params.T*params.sfreq*2,1);
   
    %get constraints
    [params.A, params.b] = getDefaultConstraints(params); 
    
    %weights
    W = [1e3; 1e3; 0; 0];
    params.W = eye(length(params.x0)).*W;
    params.Wstat = eye(length(params.x0))*1e5;
    params.R = 1;

end
