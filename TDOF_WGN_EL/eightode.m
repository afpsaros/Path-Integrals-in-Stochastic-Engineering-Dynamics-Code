function x_dot = eightode(t, x, m, c, k, e1, e2, S0, tt, BVX1, BVX2)
            
    x_dot = zeros(8, 1);

    %x1(t)=x(1)
    x_dot(1) = x(2); %x1'(t)
    x_dot(2) = x(3); %x1''(t)
    x_dot(3) = x(4); %x1'''(t)
    exp1 =...
        m^(-2)*((-3)*e1^2*k^2*x(1).^5+4*c*e2*k*x(2).^3+(-3)*c*e2*k*x(2).^2*x(...
            6)+5*c^2*x(3)+(-4)*k*m*x(3)+24*c^2*e2*x(2).^2*x(3)+15*c^2*e2^2*x(2).^...
            4*x(3)+(-6)*c^2*e2*x(2)*x(6)*x(3)+3*e1*k*x(1).^2*(k*x(5)+2*c*e2*x(2).^...
            3+c*x(6)+(-2)*m*x(3))+2*k*x(5)*(2*k+(-3)*c*e2*x(2)*x(3))+2*e1*k*x(1).^...
            3*((-4)*k+3*c*e2*x(2)*x(3))+k*x(1)*((-5)*k+(-6)*e1*m*x(2).^2+12*c*e2*...
            x(2)*x(3))+(-4)*c^2*x(7)+2*k*m*x(7)+(-3)*c^2*e2*x(2).^2*x(7));
    x_dot(4) = exp1; %x1''''(t)
    %x2(t)=x(5)
    x_dot(5) = x(6); %x2'(t)
    x_dot(6) = x(7); %x2''(t)
    x_dot(7) = x(8); %x2'''(t)
    exp2 =...
        m^(-2)*(4*k^2*x(1)+e1*k^2*x(1).^3+(-5)*k^2*x(5)+(-3)*c*e1*k*x(1).^2*x(...
            2)+c*e2*k*x(2).^3+(-4)*c^2*x(3)+2*k*m*x(3)+(-3)*c^2*e2*x(2).^2*x(3)+5*...
            c^2*x(7)+(-4)*k*m*x(7));
    x_dot(8) = exp2;
    
end