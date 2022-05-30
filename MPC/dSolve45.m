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