% A wrapper for MPC controller
%
% Parameters:
%   ps - default model & solver parameters
%   horizon - size of prediction horizon (in seconds)
%   step - size of step between predictions (in seconds)
%   Note: horizon > step and both need to divide sim time into a
%   integer of steps
function [uopt, Qstat] = MPC(ps, horizon, step)
    
    % rewrite generic params
    s_ps = ps;
    
    steps_num = ceil(ps.T/step);
    
    % unravelec control array for output
    uopt = [ps.u(1:end/2) ps.u(end/2+1:end)];

%     % solve for horizon (step 1)
%     s_ps.T = horizon;
%     % get trajectory and end conditions
%     s_ps.xref = ps.xref(1:horizon*ps.sfreq,:);
%     s_ps.xf = ps.xref(horizon*ps.sfreq,:);
%     % allocate container for control
%     s_ps.u = zeros(horizon*ps.sfreq,2);
%     [s_ps.A, s_ps.B] = getDefaultConstraints(s_ps);
    %u_insert = zeros((horizon-step)*ps.sfreq,2);
    %u_fragment = zeros((horizon)*ps.sfreq,2);

    for i = 1:steps_num
        
        fprintf("Current step: %d of %d\n", i, steps_num);

        % get simulation boundaries
        Tstart = (i-1)*step;
        Tstep = i*step;
        Tend =  min(ps.T, (i-1)*step+horizon);
        
        s_ps.T = Tend-Tstart;
        s_ps.xref = ps.xref(int32(Tstart*ps.sfreq)+1:int32(Tend*ps.sfreq)+1, :);
        s_ps.xf = ps.xref(int32(Tend*ps.sfreq)+1, :);
        s_ps.u = zeros(int32((Tend-Tstart)*ps.sfreq)*2, 1);
        [s_ps.A, s_ps.b] = getDefaultConstraints(s_ps);
      

        [u_fragment, ~] = fminconGrad(s_ps);
        
        u_sim = u_fragment(1:int32(step*ps.sfreq),:);

        % sim model and get initial position for next iter
        [~, x] = qsolve45(u_sim, Tstep-Tstart, s_ps);
        s_ps.x0 = [x(end,1:4) 0];

        uopt(int32(Tstart*ps.sfreq)+1:int32(Tstep*ps.sfreq), :) = u_sim;
       
    end
    Qstat = 0;
end

