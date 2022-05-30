function [A,b] = getDefaultConstraints(p)

    %Default constraints - u(i) < 1 && u(i) > -1
    A1 = eye(length(p.u));
    A2 = -eye(length(p.u));
    A = [A1;A2];
    b1 = ones(length(p.u),1);
    b2 = ones(length(p.u),1);
    b = [b1;b2];
end