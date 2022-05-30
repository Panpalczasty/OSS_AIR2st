function [u_qa, e_qa, pos_qa] = getQA(u, x, xref, ps)
    % uint is counted as rect approx - system is solved
    % usting step approximation of ctrl
    u_qa = 1/ps.sfreq*sum(sum(u.^2));

    % eint is counted as trapezoid approx
    err = (xref - x(:,1:4)).^2;
    e_qa = 1/ps.sfreq*(sum(sum(err)) - 0.5*(sum(err(1,:)) + sum(err(end,:))));

    % posint is counted as trapezoid approx
    pos = simpleKin(x,ps.J1, ps.J2);
    posref = simpleKin(xref, ps.J1, ps.J2);
    poserr = sqrt((pos(:,1) - posref(:,1)).^2 + (pos(:,2) - posref(:,2)).^2);
    pos_qa = 1/ps.sfreq*(sum(poserr) - 0.5*(poserr(1)+poserr(end)))/ps.T;

end