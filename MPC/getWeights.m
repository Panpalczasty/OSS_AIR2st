function ps = getWeights(ps, R, Wx, Wv, Wstat)
%UNTITLED Summary of this function goes here
% Weight manipulation
    ps.R = R;
    ps.W = [Wx;Wx;Wv;Wv].*eye(4);
    ps.Wstat = eye(4)*Wstat;

end

