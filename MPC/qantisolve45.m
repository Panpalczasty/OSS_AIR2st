%Runge_Kutty antisolver (time machine?)
function p = qantisolve45(pf, t, x, u, ps)
    %enter time machine

    %filp inputs, states and entropy increase
    u = [u; [0 0]];
    xref = flipud(ps.xref);
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
        xrh = 0.5*(xref(i+1,:) + xref(i,:));
        uh = 0.5*(u(i,:) + u(i+1,:))';

        dp1 = qantimodel(x(i,:)',    u(i+1,:)', ptmp, ps.J1, ps.J2, xref(i,:), ps.W, ps.R);    tmp = ptmp+h/2*dp1; tt = tt+h/2;
        dp2 = qantimodel(xh,         u(i+1,:)', tmp,  ps.J1, ps.J2, xrh, ps.W, ps.R);    tmp = ptmp+h/2*dp2; 
        dp3 = qantimodel(xh,         u(i+1,:)', tmp,  ps.J1, ps.J2, xrh, ps.W, ps.R);    tmp = ptmp+h*dp3;   tt = tt+h/2;
        dp4 = qantimodel(x(i+1,:)',  u(i+1,:)', tmp,  ps.J1, ps.J2, xref(i+1,:), ps.W, ps.R);

        ptmp = ptmp + h/6*(dp1 + 2*dp2 + 2*dp3 + dp4);

        %write outputs
        p(i+1,:) = ptmp';
     
    end

    %exit the time machine
    p = flipud(p);
end