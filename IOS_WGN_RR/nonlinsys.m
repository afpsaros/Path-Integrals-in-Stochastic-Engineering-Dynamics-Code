function dJdc = nonlinsys(C, g0, g1, g2, Herm, dt, m0, c0, k0, e1, e2)
                 
    y = g0 * C + Herm(:,1);
    y1 = g1 * C + Herm(:,2);
    y2 = g2 * C + Herm(:,3);

    Lhat = k0.*y+e1.*k0.*y.^3+c0.*y1+c0.*e2.*y1.^3+m0.*y2;

    % --- Calculating the integrand at discrete points ---

    dLhatdC = c0.*g1+g0.*k0+g2.*m0+3.*e1.*g0.*k0.*(Herm(:,1)+g0*C).^2+3.*c0.*e2.*g1...
                .*(Herm(:,2)+g1*C).^2;

    dJdc=(0.5 .* (Lhat(1) .* dLhatdC(1, :) + Lhat(end) .* dLhatdC(end, :)) +...
        Lhat(2:end - 1)' * dLhatdC(2:end - 1, :))' .* dt;

end
