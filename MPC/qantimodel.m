%Adjoint equations for SCARA
%addl adj variables
function dp = qantimodel(x, u, p, J1, J2, xref, W, R)

    %Steiner reduced masses
    B = J2.m*J1.d*(J2.xc*cos(x(2)) - J2.yc*sin(x(2)));
    B2 = -J2.m*J1.d*(J2.xc*sin(x(2)) + J2.yc*cos(x(2)));
    
    %inertia matrix
    M_11 = J1.I + J2.I + J2.m*J1.d^2 + 2*B;
    M_12 = B;
    M_22 = J2.I;

    %Coriolis and radial forces
    V1 = B2*(2*x(3)*x(4) + x(4)^2);
    V2 = -B2*x(3)^2;
    
    %friction 
    T1 = J1.s*tanh(10*x(3)) + J1.f*x(3);
    T2 = J2.s*tanh(10*x(4)) + J2.f*x(4);

    %build & inverse inertia matrix
    M = [M_11 M_12; M_12 M_22];
    Minv = inv(M);

    %diffs of stuff symbols : dvar
    % var - variable diffed
    %columns of derivs 

    %V1
    dV1 = [ 0;
            -B*(2*x(3)*x(4) + x(4)^2);
            +2*B2*x(4);
            +2*B2*(x(3) + x(4));
            0;0];
    
    %T1
    dT1 = [ 0;0;
            10*J1.s*1/(cosh(10*x(3))^2) + J1.f;
            0;0;0];

    %V2
    dV2 = [ 0;
            +B*x(3)^2;
            -2*B2*x(3);
            0;0;0];

    %T2
    dT2 = [ 0;0;0;
            10*J2.s*1/(cosh(10*x(4))^2) + J2.f;
            0;0];

    %inv M derivs
    dM11 = [0;
            -2*B2*(Minv(1,1)*Minv(1,2) + Minv(1,1)^2);
            0;0;0;0];
    
    dM12 = [0;
            -B2*(Minv(1,2)^2 + 2*Minv(1,2)*Minv(1,1) + Minv(1,1)*Minv(2,2))
            0;0;0;0];

    dM22 = [0;
            -2*B2*(Minv(1,2)^2 + Minv(2,2)*Minv(1,2))
            0;0;0;0];

    dH1 = [0;0;0;0;
            Minv(1,1);
            Minv(1,2)];

    dH2 = [0;0;0;0;
            Minv(1,2);
            Minv(2,2)];

    %difs of state fcn
    df1 = [0;0;1;0;0;0];
    df2 = [0;0;0;1;0;0];
    df3 = - dM11*(V1 + T1) - Minv(1,1)*(dV1 + dT1) - Minv(1,2)*(dV2 + dT2) - dM12*(V2 + T2) + dH1;
    df4 = - dM12*(V1 + T1) - Minv(1,2)*(dV1 + dT1) - Minv(2,2)*(dV2 + dT2) - dM22*(V2 + T2) + dH2;
    df5 = [0;0;0;0;0;0];
    df6 = [0;0;0;0;0;0];

    %additional dif - for cost fun SV
    dL = [W*(x(1:end-1)-xref'); +R*u];

    %difs of input fcn
    dg3 = dM11*u(1) + dM12*u(2);
    dg4 = dM12*u(1) + dM22*u(2);

    %Assemble A matrix
    A = [df1 df2 df3+dg3 df4+dg4 df5 df6];

    %condensed insanity
    %     |
    %     V
    dp =  A*p - dL;
    
    %i dont know why
    %i dont want to know why
    %but works only w/o minus

end