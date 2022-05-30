%Runge-Kutty solver
function [t,x] = qsolve45(u, tf, p)
    %rewrite inputs as horizontal
    x0 = p.x0(:);
    n = length(x0);

    %step num & size 
    Nt = round(tf*p.sfreq);
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