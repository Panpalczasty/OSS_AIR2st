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
        g(i,1) = p(i+1,5) - p(i,5);
        g(i,2) = p(i+1,6) - p(i,6);
    end
 end
