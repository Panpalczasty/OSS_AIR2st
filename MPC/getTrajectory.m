function [pos, ps] = getTrajectory(mode, ps, center, radius)
% Define trajectory here - as x = f(t), y = g(t)
% Parametric curve definition
T = ps.T;
sf = ps.sfreq;
t = linspace(0, T, T*sf);
w = 2*pi/T;

pos = zeros(length(t),2);

switch mode
    
    case 'square'
        xmin = center(1) - radius;
        xmax = center(1) + radius;
        ymin = center(2) - radius;
        ymax = center(2) + radius;
    
        pos(1:sf*T/4,1) = xmin + 8*radius*t(1:sf*T/4)/T;
        pos(1:sf*T/4,2) = ymin;
        
        pos(sf*T/4+1:sf*T/2,1) = xmax;
        pos(sf*T/4+1:sf*T/2,2) = ymin + 8*radius*t(1:sf*T/4)/T;
        
        pos(sf*T/2+1:3*sf*T/4,1) = xmax - 8*radius*t(1:sf*T/4)/T;
        pos(sf*T/2+1:3*sf*T/4,2) = ymax;
        
        pos(3*sf*T/4+1:end,1) = xmin;
        pos(3*sf*T/4+1:end,2) = ymax - 8*radius*t(1:sf*T/4)/T;

    case 'circle'
        pos(:,1) = center(1) + radius*cos(w*t);
        pos(:,2) = center(2) + radius*sin(w*t);

    case 's'
        % S shaped curve (horizontal)
        pos(1:sf*T/2,1) = radius + center(1) + radius*cos(w*t(1:sf*T/2));
        pos(1:sf*T/2,2) = center(2) + radius*sin(w*t(1:sf*T/2));
        pos(sf*T/2+1:end,1) = center(1) - radius + radius*cos(w*t(1:sf*T/2));
        pos(sf*T/2+1:end,2) = center(2) - radius*sin(w*t(1:sf*T/2));

    case 'v'
        pos(1:sf*T/2,1) = center(1) - radius + 2*radius*t(1:sf*T/2)/T;
        pos(1:sf*T/2,2) = center(2) + radius - 4*radius*t(1:sf*T/2)/T;
        pos(sf*T/2+1:end,1) = center(1) + 2*radius*t(1:sf*T/2)/T;
        pos(sf*T/2+1:end,2) = center(2) - radius + 4*radius*t(1:sf*T/2)/T;

    case 'wobble'
        pos(:,1) = center(1) + radius*cos(w*t).*(1 + 0.1*sin(10*w*t));
        pos(:,2) = center(2) + radius*sin(w*t).*(1 + 0.1*sin(10*w*t));

    case 'chirp'
        pos(:,1) = center(1) + radius*cos(10*w*t.^2/T);
        pos(:,2) = center(2) + radius*sin(10*w*t.^2/T);

    case 'agh'
        % TODO 
        pos(:,1) = 0;
        pos(:,2) = 0;
end

ps.posref = pos;
[ps.t,ps.xref] = inverseKin(pos, ps.T, "below", ps.J1, ps.J2);

ps.x0 = [ps.xref(1,1:2) 0 0 0];
ps.xf = [ps.xref(end,1:2) 0 0];
ps.xref = [ps.xref; ps.xf];

end