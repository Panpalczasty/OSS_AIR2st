
%cost function - TODO
% xf - final state 
% W = W^T >0 - weights 
function [Q,gout] = costFun(ps, u)
    
    Wstat = ps.Wstat;
    
    gout= zeros(length(u),1);
    uin = zeros(length(u)/2,2);

    %for fmincon, unravel
    uin(:,1) = u(1:end/2);
    uin(:,2) = u(end/2+1:end);

    %solve for given conditions 
    [t, x] = qsolve45(uin, ps.T, ps);
    % get state difference
    dxend = x(end,1:end-1) - ps.xf;
    % get integral - fifth state var
    x5 = x(end,end);
    % cost function
    Q = x5 + 0.5*dxend* Wstat* dxend';
    
    if nargout>1
        %get simple gradient
        g=getGrad(uin, ps);

        %for fmincon, ravel
        gout = [g(:,1); g(:,2)];
    end
end