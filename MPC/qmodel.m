%SCARA 2axis robot model
%with additional state eq - for cost function
%dx5 = (x - xref)'*W*(x - xref)
function dx = qmodel(x, u, J1, J2, xref, W, R)
    x = x(:);
    %calculate additional params 
    
    %Steiner reduced masses
    B = J2.m*J1.d*(J2.xc*cos(x(2)) - J2.yc*sin(x(2)));
    B2 = -J2.m*J1.d*(J2.xc*sin(x(2)) + J2.yc*cos(x(2)));
    
    %inertia matrix
    M_11 = J1.I + J2.I + J2.m*J1.d^2+2*B;
    M_12 = B;
    M_22 = J2.I;

    %Coriolis and radial forces
    V1 = B2*(2*x(3)*x(4)+ x(4)^2);
    V2 = -B2*x(3)^2;
    
    %friction 
    T1 = J1.s*tanh(10*x(3)) + J1.f*x(3);
    T2 = J2.s*tanh(10*x(4)) + J2.f*x(4);

    %build & inverse inertia matrix
    M = [M_11 M_12; M_12 M_22];
    Minv = inv(M);
    
    %build equations
    f0 = [  x(3); 
            x(4); 
            -Minv(1,1)*(V1+T1) - Minv(1,2)*(V2+T2);
            -Minv(1,2)*(V1+T1) - Minv(2,2)*(V2+T2);
            0.5*(x(1:end-1)' - xref)*W*(x(1:end-1) - xref') + 0.5*R*u(1)^2 + 0.5*R*u(2)^2
         ];

    g1 = [  0;
            0;
            Minv(1,1);
            Minv(1,2);
            0
         ];
    
    g2 = [  0;
            0;
            Minv(1,2);
            Minv(2,2);
            0
         ];
    
    %calculate output
    dx = f0 + g1*u(1) + g2*u(2);
end
