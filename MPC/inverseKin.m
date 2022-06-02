%inverse kinematics - above or below
function [t, x] = inverseKin(pos, T, config, J1, J2)

    rad = sqrt(pos(:,1).^2 + pos(:,2).^2);
    b = (rad.^2 - J1.d^2 - J2.d^2)/(2*J1.d*J2.d);

    samples = length(pos(:,1));
    sfreq = samples(1)/T;
    t = linspace(0, T, samples);

    if config == "above"
       c = 1; 
    elseif config == "below"
       c = -1;
    end
 
    if any(abs(b)>1)  
        disp("Trajectory out of robot range!");
        x = zeros(samples,4);
    else
        alpha = atan2(-c*J2.d*sqrt(1-b.^2), J1.d + J2.d*b);
    
        x1 = alpha + atan2(pos(:,2), pos(:,1));
        x2 = c*acos(b); 
        x3 = [diff(x1)*sfreq; 0];
        x4 = [diff(x2)*sfreq; 0];    
    
        x3(end) = x3(end-1);
        x4(end) = x4(end-1);

        x = [x1 x2 x3 x4];
    end
end
