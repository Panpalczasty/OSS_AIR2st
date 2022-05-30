%get vector version of Hu
function Hu = getHuvec(t,x,p,u,params)
    Hu = zeros(length(t),2);
    u = [u; u(end,:)];
    for i = 1:length(t)
        Hu(i,:) = getHu(t, p(i,:), x(i,:), u(i,:), params.J1, params.J2, params.R);
    end
end