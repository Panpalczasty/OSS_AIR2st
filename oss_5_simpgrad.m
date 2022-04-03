clc
clear all
close all

%% Desc

%This file implements optimization for changes in input value
%at a discrete raster - with zoh discrete changes
%input changest with the same sampling freq as sim
%cost fun is a combination of energy minimalization, trajectory deviation
%and final position deviation

%% get model parameters
%no switch times, only levels here

params = getDefaultParams();

%Weight manipulation
W = 1*[3e1; 3e1; 1; 1];
Wstat = 1*[1e3; 1e3; 1e2; 1e2];

params.W = eye(length(params.x0)).*W;
params.Wstat = eye(length(params.x0)).*Wstat;
params.R = 0;

params.sfreq = 20;
params.T = 2.5;
T = params.T;
%Define trajectory here - as x = f(t), y = g(t)
t = linspace(0, params.T, params.T*params.sfreq);
w = 2*pi/T;
pos(:,1) = 0.1*cos(0.25*w*t)+1;
pos(:,2) = 0.1*sin(0.25*w*t); 
params.posref = pos;

[params.t,params.xref] = inverseKin(pos, params.T, "below", params.J1, params.J2);
params.x0 = [params.xref(1,:) 0];
params.xf = params.xref(end,:);
params.u = zeros((length(t)),2);


%% optimizers - pick one

alpha = 1e-3;
iter = 2500;

[uopt, Qstat] = simpleGrad(params, alpha, iter);
%[uopt, Qstat] = iterGrad(params, alpha, iter);
%[uopt, Qstat] = fminconGrad(params);

%% get stats

disp(costFun(params,uopt))
 
%solve
[t, x] = qsolve45(uopt, params.T, params);
%antisolve
pf =  - (params.Wstat)*(x(end,1:end-1) - params.xf)';
pf = [pf;0;0];
p = qantisolve45(pf, t, x, uopt, params);

pos = simpleKin(x, params.J1, params.J2);
Hu = getHuvec(t, x, p, uopt, params);

%% perform checks 

%numerical gradient = analitycal gradient
%checkGradient(params, 0.01);

% psi(0) = - derivative of cost fun w.r.t state variables
%checkAdjoint(p, params, 0.01);

%% plots, plots, plots
plotAdjoint(1,t,p,params);
plotAngles(2,t,x,uopt,params);
plotPos(3, pos, params.posref, params);
plotSwitch(5,t,uopt,Hu,params);
plotCost(4,Qstat,params);

animate(6,pos, params);

%% functions - core

%SCARA 2axis robot model
%with additional state eq - for cost function
%dx5 = (x - xref)'*Q*(x - xref)
function dx = qmodel(x, u, J1, J2, xref, W, R)
    x = x(:);
    %calculate additional params 
    
    %Steiner reduced masses
    B = J2.m*J1.d*(J2.xc*cos(x(2)) - J2.yc*sin(x(2)));
    B2 = -J2.m*J1.d*(J2.xc*sin(x(2)) + J2.yc*cos(x(2)));
    
    %inertia matrix
    M_11 = J1.I + J2.I + J2.m*J1.d^2+2*B;
    M_12 = B;
    M_22 = J2.I;

    %Coriolis and radial forces
    V1 = B2*(2*x(3)*x(4)+ x(4)^2);
    V2 = -B2*x(3)^2;
    
    %friction 
    T1 = J1.s*tanh(10*x(3)) + J1.f*x(3);
    T2 = J2.s*tanh(10*x(4)) + J2.f*x(4);

    %build & inverse inertia matrix
    M = [M_11 M_12; M_12 M_22];
    Minv = inv(M);
    
    %build equations
    f0 = [  x(3); 
            x(4); 
            -Minv(1,1)*(V1+T1) - Minv(1,2)*(V2+T2);
            -Minv(1,2)*(V1+T1) - Minv(2,2)*(V2+T2);
            0.5*(x(1:end-1)' - xref)*W*(x(1:end-1) - xref') + 0.5*R*u(1)^2 + 0.5*R*u(2)^2
         ];

    g1 = [  0;
            0;
            Minv(1,1);
            Minv(1,2);
            0
         ];
    
    g2 = [  0;
            0;
            Minv(1,2);
            Minv(2,2);
            0
         ];
    
    %calculate output
    dx = f0 + g1*u(1) + g2*u(2);
end

%Adjoint equations for SCARA
%addl adj variable
function dp = qantimodel(x, u, p, J1, J2, xref, W, R)

    %Steiner reduced masses
    B = J2.m*J1.d*(J2.xc*cos(x(2)) - J2.yc*sin(x(2)));
    B2 = -J2.m*J1.d*(J2.xc*sin(x(2)) + J2.yc*cos(x(2)));
    
    %inertia matrix
    M_11 = J1.I + J2.I + J2.m*J1.d^2 + 2*B;
    M_12 = B;
    M_22 = J2.I;

    %Coriolis and radial forces
    V1 = B2*(2*x(3)*x(4) + x(4)^2);
    V2 = -B2*x(3)^2;
    
    %friction 
    T1 = J1.s*tanh(10*x(3)) + J1.f*x(3);
    T2 = J2.s*tanh(10*x(4)) + J2.f*x(4);

    %build & inverse inertia matrix
    M = [M_11 M_12; M_12 M_22];
    Minv = inv(M);

    %diffs of stuff symbols : dvar
    % var - variable diffed
    %columns of derivs 

    %V1
    dV1 = [ 0;
            -B*(2*x(3)*x(4) + x(4)^2);
            +2*B2*x(4);
            +2*B2*(x(3) + x(4));
            0;0];
    
    %T1
    dT1 = [ 0;0;
            10*J1.s*1/(cosh(10*x(3))^2) + J1.f;
            0;0;0];

    %V2
    dV2 = [ 0;
            +B*x(3)^2;
            -2*B2*x(3);
            0;0;0];

    %T2
    dT2 = [ 0;0;0;
            10*J2.s*1/(cosh(10*x(4))^2) + J2.f;
            0;0];

    %inv M derivs
    dM11 = [0;
            -2*B2*(Minv(1,1)*Minv(1,2) + Minv(1,1)^2);
            0;0;0;0];
    
    dM12 = [0;
            -B2*(Minv(1,2)^2 + 2*Minv(1,2)*Minv(1,1) + Minv(1,1)*Minv(2,2))
            0;0;0;0];

    dM22 = [0;
            -2*B2*(Minv(1,2)^2 + Minv(2,2)*Minv(1,2))
            0;0;0;0];

    dH1 = [0;0;0;0;
            Minv(1,1);
            Minv(1,2)];

    dH2 = [0;0;0;0;
            Minv(1,2);
            Minv(2,2)];

    %difs of state fcn
    df1 = [0;0;1;0;0;0];
    df2 = [0;0;0;1;0;0];
    df3 = - dM11*(V1 + T1) - Minv(1,1)*(dV1 + dT1) - Minv(1,2)*(dV2 + dT2) - dM12*(V2 + T2) + dH1;
    df4 = - dM12*(V1 + T1) - Minv(1,2)*(dV1 + dT1) - Minv(2,2)*(dV2 + dT2) - dM22*(V2 + T2) + dH2;
    df5 = [0;0;0;0;0;0];
    df6 = [0;0;0;0;0;0];

    %additional dif - for cost fun SV
    dL = [W*(x(1:end-1)-xref'); +R*u];

    %difs of input fcn
    dg3 = dM11*u(1) + dM12*u(2);
    dg4 = dM12*u(1) + dM22*u(2);

    %Assemble A matrix
    A = [df1 df2 df3+dg3 df4+dg4 df5 df6];

    %condensed insanity
    %     |
    %     V
    dp =  A*p - dL;
    
    %i dont know why
    %i dont want to know why
    %but works only w/o minus

end

%Runge_Kutty antisolver (time machine?)
function p = qantisolve45(pf, t, x, u, ps)
    %enter time machine

    %filp inputs, states and entropy increase
    u = [u; [0 0]];
    x = flipud(x);
    u = flipud(u);

    %reverse temporal flow
    t = t(end) - t;
    t = flipud(t);

    %configure spacetime sampling
    Nt = length(t);
    h = tf/Nt;

    %declare degrees of freedom 
    n = length(pf);
    p = zeros(Nt, n);
    p(1,:) = pf';
   
    %finalize preparations
    ptmp = pf;
    tt = 0;

    tmp = zeros(n,1);
    dp1 = zeros(n,1);
    dp2 = zeros(n,1);
    dp3 = zeros(n,1);
    dp4 = zeros(n,1);
    
    %perform time travel
    for i = 1:Nt-1
        h = t(i+1) - t(i);
        xh = 0.5*(x(i,:) + x(i+1,:))';

        dp1 = qantimodel(x(i,:)',    u(i+1,:)', ptmp, ps.J1, ps.J2, ps.xref(i,:), ps.W, ps.R);    tmp = ptmp+h/2*dp1; tt = tt+h/2;
        dp2 = qantimodel(xh,         u(i+1,:)', tmp,  ps.J1, ps.J2, ps.xref(i,:), ps.W, ps.R);    tmp = ptmp+h/2*dp2; 
        dp3 = qantimodel(xh,         u(i+1,:)', tmp,  ps.J1, ps.J2, ps.xref(i,:), ps.W, ps.R);    tmp = ptmp+h*dp3;   tt = tt+h/2;
        dp4 = qantimodel(x(i+1,:)',  u(i+1,:)', tmp,  ps.J1, ps.J2, ps.xref(i,:), ps.W, ps.R);

        ptmp = ptmp + h/6*(dp1 + 2*dp2 + 2*dp3 + dp4);

        %write outputs
        p(i+1,:) = ptmp';
     
    end

    %exit the time machine
    p = flipud(p);
end

%Runge-Kutty solver
function [t,x] = qsolve45(u, tf, p)
    %rewrite inputs as horizontal
    x0 = p.x0(:);
    n = length(x0);

    %step num & size 
    Nt = ceil(tf*p.sfreq);
    h = tf/Nt;

    %declare matrices
    x = zeros(Nt+1, n);
    t = zeros(Nt+1, 1);
    x(1,:) = x0';

    xtmp = x0;
    tt = 0;

    tmp = zeros(n,1);
    dx1 = zeros(n,1);
    dx2 = zeros(n,1);
    dx3 = zeros(n,1);
    dx4 = zeros(n,1);
    
    %solve 
    for i = 1:Nt
        dx1 = qmodel(xtmp, u(i,:), p.J1, p.J2, p.xref(i,:), p.W, p.R);   tmp = xtmp+h/2*dx1; tt = tt+h/2;
        dx2 = qmodel(tmp, u(i,:), p.J1, p.J2, p.xref(i,:), p.W, p.R);    tmp = xtmp+h/2*dx2; 
        dx3 = qmodel(tmp, u(i,:), p.J1, p.J2, p.xref(i,:), p.W, p.R);    tmp = xtmp+h*dx3;   tt = tt+h/2;
        dx4 = qmodel(tmp, u(i,:), p.J1, p.J2, p.xref(i,:), p.W, p.R);

        xtmp = xtmp + h/6*(dx1 + 2*dx2 + 2*dx3 + dx4);

        %write outputs
        x(i+1,:) = xtmp';
        t(i+1) = tt;
    end
end

%ode45 solver for discontinuities -TODO
%i dont really need to do this
function [t,x,uk,nseg] = dSolve45(p, u)

    [ntau, nu] = size(u);
    n = length(p.x0);         
    tau = [0; tau];

    %Get num of samples
    Nt = 1+sum(ceil(p.sfreq*diff(tau)));

    %Make output vector
    x = zeros(Nt, n);    
    uk = zeros(Nt, nu);
    t = zeros(Nt, 1);
    nseg = zeros(ntau,1);
    k = 1;                  %starting sample

    for i = 1:ntau-1 %for i-th continuous segment
        simt = tau(i+1) - tau(i); %get partial sim time
        Ni = ceil(p.sfreq*simt); %set samples num for current sim
        ui = kron(u(i,:),ones(Ni,1));
        [tseg,xseg] = solve45(ui, simt, p); %solve with i-th input
        p.x0 = xseg(end,:)'; %get next boundary
        x(k:k+Ni-1,:) = xseg(1:end-1,:); %write new states
        t(k:k+Ni-1) = tau(i) + tseg(1:end-1,:); %write new time
        uk(k:k+Ni-1,:) = ui;
        k = k+Ni; 
        nseg(i) = k;
    end

    x(end,:) = xseg(end,:); %set end state
    uk(end,:) = u(end,:); %set end input
    t(end) = tau(end); %set end time

end

%cost function - TODO
% xf - final state 
% W = W^T >0 - weights 
function [Q,g] = costFun(ps, u)
    
    Wstat = ps.Wstat;

    %solve for given conditions 
    [t, x] = qsolve45(u, ps.T, ps);
    % get state difference
    dxend = x(end,1:end-1) - ps.xf;
    % get integral - fifth state var
    x5 = x(end,end);
    % cost function
    Q = x5 + 0.5*dxend* Wstat* dxend';
    
    if nargout>1
        %get simple gradient
        g=getGrad(u, ps);
    end
end

%get hamiltionian gradients
function g = getGrad(u, ps)
    nu = length(u(:,1));
    %go forward in time and observe model
    [t, x] = qsolve45(u, ps.T, ps);

    %get time machine
    pf = - ps.Wstat*(x(end,1:end-1) - ps.xf)';
    pf = [pf;0;0];
        
    %go back in time and de-observe model
    p = qantisolve45(pf, t, x, u, ps);
    
    %get your gradient container primed and ready
    g = zeros(nu, 2);
    
    %and calculate it like there's no tomorrow
    for i = 1: nu
        %hamiltionian simple gradient
        %g(i,:) = getHu(t,p,x,u,ps.J1,ps.J2,ps.R)/ps.sfreq;
        g(i,1) =  p(i+1,5) - p(i,5);
        g(i,2) =  p(i+1,6) - p(i,6);
    end
 end

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

%get switch fcn
function Hu = getHu(t, p, x, u, J1, J2, R)
    
    %Steiner reduced masses
    B = J2.m*J1.d*(J2.xc*cos(x(2)) - J2.yc*sin(x(2)));
    B2 = -J2.m*J1.d*(J2.xc*sin(x(2)) + J2.yc*cos(x(2)));
    
    %inertia matrix
    M_11 = J1.I + J2.I + J2.m*J1.d^2+2*B;
    M_12 = B;
    M_22 = J2.I;
    
    %build & inverse inertia matrix
    M = [M_11 M_12; M_12 M_22];
    Minv = inv(M);
    
    Hu = [Minv(1,1)*p(3) + Minv(1,2)*p(4) - R*u(1)
          Minv(1,2)*p(3) + Minv(2,2)*p(4) - R*u(2)];
end


%% functions - utilities

%declare modelfunction [t, x] = inverseKin(pos, T, config, J1, J2)
function params = getDefaultParams()

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
    
    %Default sampling frequency
    params.sfreq = 100;
    params.T = 10;

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

    %get constraints
    %[params.A, params.b] = getDefaultConstraints(params); 

    %initial input - const
    params.u = zeros(params.T*params.sfreq,2);
    
    %weights
    W = [1e3; 1e3; 0; 0];
    params.W = eye(length(params.x0)).*W;
    params.Wstat = eye(length(params.x0))*1e5;
    params.R = 1;

end

function [A,b] = getDefaultConstraints(p)

    %Default
    nb = length(p.u(:,1));
    b = 0.0*ones(nb+1,1);

end

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

%simple kinematics - visualize robot pos
function pos = simpleKin(Jdata, J1, J2)
    N = length(Jdata(:,1));
    pos = zeros(N,4);
    for i = 1:N
        Rx1 = [cos(Jdata(i,1)) -sin(Jdata(i,1));
               sin(Jdata(i,1)) cos(Jdata(i,1))];
        Rx2 = [cos(Jdata(i,2)) -sin(Jdata(i,2));
               sin(Jdata(i,2)) cos(Jdata(i,2))];
        pos(i,1:2) = Rx1*[J1.d;0] + Rx1*Rx2*[J2.d; 0];
        pos(i,3:4) = Rx1*[J1.d;0];
    end
end

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
    
        x = [x1 x2 x3 x4];
    end
end

%gradient based on fmincon - TODO
function [uopt,Qstat] = fminconGrad(ps)

    Q=@(u) costFun(ps,u);
    
    options = optimoptions('fmincon');
    options.SpecifyObjectiveGradient = true;
    options.Display = 'iter';
    options.StepTolerance = 1e-5;
    options.Algorithm = 'interior-point';

    uopt = fmincon(Q,ps.u,[],[],[],[],[],[],[],options);
    
    Qstat = costFun(ps, uopt);
end

%get vector version of Hu
function Hu = getHuvec(t,x,p,u,params)
    Hu = zeros(length(t),2);
    u = [u; u(end,:)];
    for i = 1:length(t)
        Hu(i,:) = getHu(t, p(i,:), x(i,:), u(i,:), params.J1, params.J2, params.R);
    end
end

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


%% functions - visualization 

%animate robot
function res = animate(id, pos, ps)

    N = length(pos(:,1));
    figure(id);
    grid();
    hold on;
    plot(0,0,"k*",'LineWidth',3);
    j2plot = plot([pos(1,3) pos(1,1)],[pos(1,4) pos(1,2)], 'k*-', 'LineWidth',3);
    j1plot = plot([0 pos(1,3)],[0 pos(1,4)],'k','LineWidth',3);
    axis([-2 2 -2 2]);
    
    xlabel("X [m]");
    ylabel("Y [m]");
    title("Robot Animation")
    
    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);

    filename = sprintf("gifs/%.2f_%.2f_to_%.2f_%.2f.gif", x01, x02, xf1, xf2);

    for i = 1:1:N
        j1plot.XData = [0 pos(i,3)];
        j1plot.YData = [0 pos(i,4)];
        j2plot.XData = [pos(i,3) pos(i,1)];
        j2plot.YData = [pos(i,4) pos(i,2)];
        drawnow 
        
        frame = getframe(id);
        img = frame2im(frame);
        [img, cm] = rgb2ind(img,256);

        if i == 1
            imwrite(img, cm, filename, 'gif', 'LoopCount',inf);
        else
            imwrite(img, cm, filename, 'gif', "WriteMode","append", 'DelayTime', 1/(2*ps.sfreq));
        end
    end
    
    res = 0;
end

%plot angles in time domain
function res = plotAngles(id, tout, xout, uout, ps)

    f = figure(id);
    f.Position = [600 0 1000 900];
    uout = [uout; uout(end,:)];

    %plot joint values
    subplot(3,1,1);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Angular position [rad]");
    title("Position of J1, J2");
    plot(tout,xout(:,1),"b-","LineWidth",1.5);
    plot(tout,xout(:,2),"r-","LineWidth",1.5);
    legend("$x_1$","$x_2$","Interpreter","latex");
    yticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2 3*pi/4 pi]);
    yticklabels({'-\pi','-0.75\pi','-0.5\pi','-0.25\pi','0','0.25\pi','0.5\pi','0.75\pi','\pi'});

    %plot joint velocities
    subplot(3,1,2);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Angular velocity [rad/s]");
    title("Angular velocity of J1,J2");
    plot(tout,xout(:,3),"b-","LineWidth",1.5);
    plot(tout,xout(:,4),"r-","LineWidth",1.5);
    legend("$\dot{x_1}$", "$\dot{x_2}$","Interpreter","latex")

    %plot input
    subplot(3,1,3);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Input Value");
    title("Input");
    plot(tout,uout(:,1),"b-","LineWidth",1.5);
    plot(tout,uout(:,2),"r-","LineWidth",1.5);
    %ylim([-1.2 1.2])
    legend("$u_1$","$u_2$","Interpreter","latex")
    
    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);
    filename = sprintf("plots/state/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    saveas(id,filename);

    res = 0;
end

%plot time travel parameters
function res = plotAdjoint(id, tout, pout, ps)
    
    f = figure(id);
    f.Position = [600 0 1000 600];

    %plot position adjoints
    subplot(2,1,1);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Weight component");
    title("Adjoint of position of J1, J2");
    plot(tout,pout(:,1),"b--","LineWidth",1.5);
    plot(tout,pout(:,2),"r--","LineWidth",1.5);
    legend("$\psi_1$","$\psi_2$","Interpreter","latex")

    %plot joint velocities
    subplot(2,1,2);
    hold on;
    grid on;
    xlabel("Time [s]");
    ylabel("Weight component");
    title("Adjoint of velocity of J1,J2");
    plot(tout,pout(:,3),"b--","LineWidth",1.5);
    plot(tout,pout(:,4),"r--","LineWidth",1.5);
    legend("$\psi_3$", "$\psi_4$", "Interpreter","latex")

    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);
    filename = sprintf("plots/adj/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    saveas(id,filename);

    res = 0;
end

%plot trajectory in 2D plane
function res = plotPos(id, pos, posref, ps)
    
    f = figure(id);
    f.Position = [0 600 600 400];


    plot(pos(:,1),pos(:,2), "k","LineWidth",1.2);
    hold on;
    plot(posref(:,1), posref(:,2),"b--");
    grid();
    plot(pos(1,1),pos(2,1), "b*");
    xlabel("x[m]");
    ylabel("y[m]");
    title("Trajectory of robot");
    legend("Trajectory", "Reference trajectory", "Start point");

    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);
    filename = sprintf("plots/pos/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    %saveas(id,filename);

    res = 0;
end

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

%plot switch
function res = plotSwitch(id, tout, uout, Huout, ps)
    
    uout = [uout; uout(end,:)];
    f = figure(id);
    f.Position = [0 0 600 400];

    subplot(2,1,1);
    plot(tout, uout(:,1),"b-", "LineWidth", 1.5);
    hold on;
    plot(tout, Huout(:, 1), "b--");

    grid;
    xlabel("Time [s]");
    ylabel("Value");
    title("Switch fcn - JT1");
    labels = ["$u_1$", "$\phi_1$"];
    legend(labels, "Interpreter","latex");

    subplot(2,1,2);
    plot(tout, uout(:,2), "r-", "LineWidth", 1.5);
    hold on;
    plot(tout, Huout(:, 2), "r--");

    grid;
    xlabel("Time [s]");
    ylabel("Value");
    title("Switch fcn - JT2");
    labels = ["$u_2$", "$\phi_2$"];
    legend(labels, "Interpreter","latex");

    x01 = ps.x0(1); 
    x02 = ps.x0(2);
    xf1 = ps.xf(1);
    xf2 = ps.xf(2);
    filename = sprintf("plots/switch/%.2f_%.2f_to_%.2f_%.2f.png", x01, x02, xf1, xf2);
    saveas(id,filename);
 
    res = 0;
end


