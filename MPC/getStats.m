function [x,t,p,hu,pos] = getStats(ps, uopt)

    % solve with optimized control 
    [t, x] = qsolve45(uopt, ps.T, ps);
    
    % antisolve with optimized control
    pf =  - (ps.Wstat)*(x(end,1:end-1) - ps.xf)';
    pf = [pf;0;0];
    p = qantisolve45(pf, t, x, uopt, ps);
    
    % get position in cartesian coords
    pos = simpleKin(x, ps.J1, ps.J2);
    
    % get hamiltonian deriv wrt control
    hu = getHuvec(t, x, p, uopt, ps);

end