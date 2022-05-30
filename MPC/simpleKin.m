%simple kinematics - visualize robot pos
function pos = simpleKin(Jdata, J1, J2)
    N = length(Jdata(:,1));
    pos = zeros(N,4);
    for i = 1:N
        Rx1 = [cos(Jdata(i,1)) -sin(Jdata(i,1));
               sin(Jdata(i,1)) cos(Jdata(i,1))];
        Rx2 = [cos(Jdata(i,2)) -sin(Jdata(i,2));
               sin(Jdata(i,2)) cos(Jdata(i,2))];
        pos(i,1:2) = Rx1*[J1.d;0] + Rx1*Rx2*[J2.d; 0];
        pos(i,3:4) = Rx1*[J1.d;0];
    end
end