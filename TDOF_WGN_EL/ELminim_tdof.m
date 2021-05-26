function val = ELminim_tdof(m, c, k, e1, e2, t, N, grid, S0, numt, numinit, ...
        maxiter, tol, setopt, parallel)   
    
    val = NaN(N, 1);
    if parallel == 1
        parfor_progress(N);
        parfor r = 1:N
            tf = t(r);    

            BVX1 = [grid(r, 1) grid(r, 3)];
            BVX2 = [grid(r, 2) grid(r, 4)];
            
            tt = linspace(0, tf, numt); % time-domain discretization
              
            if setopt == 1
                options = bvpset('NMax', maxiter, 'RelTol', tol);
            else
                options = []; 
            end
            
            solinit = bvpinit(linspace(0, tt(end), numinit), zeros(1, 8));
            sol = bvp4c(@eightode, @eightbc, solinit, options, m, c, k, e1, e2, S0, tt, BVX1, BVX2);
            XX = deval(sol, tt); 
            Integ = zeros(1, length(tt));
            Integ = ...
                (1/4)*pi^(-1)*S0^(-1)*((2*k*XX(1,:)+e1*k*XX(1,:).^3+(-1)*k*XX(5,:)+2*...
                    c*XX(2,:)+c*e2*XX(2,:).^3+(-1)*c*XX(6,:)+m*XX(3,:)).^2+(k*XX(1,:)+(-2)*...
                    k*XX(5,:)+c*XX(2,:)+(-2)*c*XX(6,:)+(-1)*m*XX(7,:)).^2); 
            
            result = trapz(tt, Integ); 
            val(r) = result;

            parfor_progress;
        end
        parfor_progress(0);
    else
        for r = 1:N
            
            r
            tf = t(r);    

            BVX1 = [grid(r, 1) grid(r, 3)];
            BVX2 = [grid(r, 2) grid(r, 4)];
            
            tt = linspace(0, tf, numt); % time-domain discretization
              
            if setopt == 1
                options = bvpset('NMax', maxiter, 'RelTol', tol);
            else
                options = []; 
            end
            
            solinit = bvpinit(linspace(0, tt(end), numinit), zeros(1, 8));
            sol = bvp4c(@eightode, @eightbc, solinit, options, m, c, k, e1, e2, S0, tt, BVX1, BVX2);
            XX = deval(sol, tt); 
            Integ = zeros(1, length(tt));
            Integ = ...
                (1/4)*pi^(-1)*S0^(-1)*((2*k*XX(1,:)+e1*k*XX(1,:).^3+(-1)*k*XX(5,:)+2*...
                    c*XX(2,:)+c*e2*XX(2,:).^3+(-1)*c*XX(6,:)+m*XX(3,:)).^2+(k*XX(1,:)+(-2)*...
                    k*XX(5,:)+c*XX(2,:)+(-2)*c*XX(6,:)+(-1)*m*XX(7,:)).^2); 
          
            result = trapz(tt, Integ); 
            val(r) = result;

        end   
    end
end