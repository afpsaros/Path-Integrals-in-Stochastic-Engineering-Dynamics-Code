function [action, elapsed_time] = ELminim(system_params, eval_params, algo_params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m0 = system_params.m0;
    c0 = system_params.c0;
    k0 = system_params.k0;
    S0 = system_params.S0;
    e1 = system_params.e1;
    e2 = system_params.e2;
    p0 = system_params.p0;
    q0 = system_params.q0;
    r0 = system_params.r0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ti = eval_params.ti;
    grid = eval_params.grid;
    N = eval_params.N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    numt = algo_params.numt;
    numinit = algo_params.numinit;
    maxiter = algo_params.maxiter;
    tol = algo_params.tol;
    parallel = algo_params.parallel;
    setopt = algo_params.setopt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    action = NaN(N, 1);
    tic
    if parallel == 1
        parfor_progress(N);
        parfor r = 1:N
            tf = grid(r, 5);    

            BVX = grid(r, 1:4);
            
            tt = linspace(ti, tf, numt); % time-domain discretization
              
            if setopt == 1
                options = bvpset('NMax', maxiter, 'RelTol', tol);
            else
                options = []; 
            end
            
            solinit = bvpinit(linspace(ti, tt(end), numinit), zeros(1, 8));
            sol = bvp4c(@eightode, @eightbc, solinit, options, m0, c0, k0, e1, e2, p0, q0, r0, S0, ...
                tt, BVX);
            XX = deval(sol, tt); 
            lagrangian = zeros(1, length(tt));
            lagrangian = ...
                    (1/4).*pi^(-1).*S0^(-1).*(r0.*(k0.*XX(1,:)+e1.*k0.*XX(1,:).^3+c0.*XX(...
                    2,:)+c0.*e2.*XX(2,:).^3+m0.*XX(3,:))+q0.*(k0.*XX(2,:)+3.*e1.*k0.*XX(1,...
                    :).^2.*XX(2,:)+c0.*XX(3,:)+3.*c0.*e2.*XX(2,:).^2.*XX(3,:)+m0.*XX(4,:))...
                    +p0.*(k0.*XX(3,:)+e1.*k0.*(6.*XX(1,:).*XX(2,:).^2+3.*XX(1,:).^2.*XX(3,...
                    :))+c0.*XX(4,:)+c0.*e2.*(6.*XX(2,:).*XX(3,:).^2+3.*XX(2,:).^2.*XX(4,:)...
                    )+m0.*XX(5,:))).^2;
          
            action(r) = trapz(tt, lagrangian); 

            parfor_progress;
        end
        parfor_progress(0);
    else
        for r = 1:N
            
            r
            tf = grid(r, 5);    

            BVX = grid(r, 1:4);
            
            tt = linspace(ti, tf, numt); % time-domain discretization
              
            if setopt == 1
                options = bvpset('NMax', maxiter, 'RelTol', tol);
            else
                options = []; 
            end
            
            solinit = bvpinit(linspace(ti, tt(end), numinit), zeros(1, 8));
            sol = bvp4c(@eightode, @eightbc, solinit, options, m0, c0, k0, e1, e2, p0, q0, r0, S0, ...
                tt, BVX);
            XX = deval(sol, tt); 
            lagrangian = zeros(1, length(tt));
            lagrangian = ...
                    (1/4).*pi^(-1).*S0^(-1).*(r0.*(k0.*XX(1,:)+e1.*k0.*XX(1,:).^3+c0.*XX(...
                    2,:)+c0.*e2.*XX(2,:).^3+m0.*XX(3,:))+q0.*(k0.*XX(2,:)+3.*e1.*k0.*XX(1,...
                    :).^2.*XX(2,:)+c0.*XX(3,:)+3.*c0.*e2.*XX(2,:).^2.*XX(3,:)+m0.*XX(4,:))...
                    +p0.*(k0.*XX(3,:)+e1.*k0.*(6.*XX(1,:).*XX(2,:).^2+3.*XX(1,:).^2.*XX(3,...
                    :))+c0.*XX(4,:)+c0.*e2.*(6.*XX(2,:).*XX(3,:).^2+3.*XX(2,:).^2.*XX(4,:)...
                    )+m0.*XX(5,:))).^2;
          
            action(r) = trapz(tt, lagrangian); 

        end   
    end
    elapsed_time = toc;
end