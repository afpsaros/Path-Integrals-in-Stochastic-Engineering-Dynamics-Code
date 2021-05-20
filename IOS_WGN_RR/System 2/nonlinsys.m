function dJdc = nonlinsys(C, g0, g1, g2, Herm, dt, m0, c0, k0, e1, lambda)
                 
    y = g0 * C + Herm(:,1);
    y1 = g1 * C + Herm(:,2);
    y2 = g2 * C + Herm(:,3);

    Lhat = k0.*y+c0.*y1+m0.*y2+e1.*y1.*abs(y1);

    % --- Calculating the integrand at discrete points ---

    dLhatdC = c0.*g1+g0.*k0+g2.*m0+e1.*g1.*abs(Herm(:,2)+g1*C)+2.*e1.*g1.*pi.^(-1).*(Herm(:,2)+g1*C).*atan(lambda.*(Herm(:,2)+g1*C));

    dJdc=(0.5 .* (Lhat(1) .* dLhatdC(1, :) + Lhat(end) .* dLhatdC(end, :)) +...
        Lhat(2:end - 1)' * dLhatdC(2:end - 1, :))' .* dt;

end
