%get switch fcn
function Hu = getHu(t, p, x, u, J1, J2, R)
    
    %Steiner reduced masses
    B = J2.m*J1.d*(J2.xc*cos(x(2)) - J2.yc*sin(x(2)));
    B2 = -J2.m*J1.d*(J2.xc*sin(x(2)) + J2.yc*cos(x(2)));
    
    %inertia matrix
    M_11 = J1.I + J2.I + J2.m*J1.d^2+2*B;
    M_12 = B;
    M_22 = J2.I;
    
    %build & inverse inertia matrix
    M = [M_11 M_12; M_12 M_22];
    Minv = inv(M);
    
    Hu = [Minv(1,1)*p(3) + Minv(1,2)*p(4) - R*u(1)
          Minv(1,2)*p(3) + Minv(2,2)*p(4) - R*u(2)];
end