function [action, exit_flag, elapsed_time] = RRminim(system_params, eval_params, algo_params)   

    action = [];
    exit_flag = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m0 = system_params.m0;
    c0 = system_params.c0;
    k0 = system_params.k0;
    S0 = system_params.S0;
    e1 = system_params.e1;
    lambda = system_params.lambda;
    ord = system_params.ord;
    ndof = system_params.ndof;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ti = eval_params.ti;
    grid = eval_params.grid;
    N = eval_params.N;
    Nt = eval_params.Nt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    h = algo_params.h;
    Fs = algo_params.Fs;
    parallel = algo_params.parallel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%
    % Boundary conditions
    con = zeros(ndof, 4, N/Nt); % 4 because we have 2 initial constraints and 2 final
    con(1,3,:) = grid(1:N/Nt, 1);
    con(1,4,:) = grid(1:N/Nt, 2);
    
    tic
    for jj = 1:Nt

        tf = grid(jj + (jj - 1) * N / Nt, 3);

        % - Construction shifted Legendre polynomials -
        syms t
        P = sym('P', [1 h]);
        P(1) = 1;
        P(2) = (2 * t - ti - tf) / (tf - ti);

        for kk = 3:h
            p = kk - 2;
            P(kk) = expand((2 * p + 1) / (p + 1) * (2 * t - ti - tf)...
                / (tf - ti) * P(kk - 1) - p / (p + 1) * P(kk - 2));
        end
        clear kk p

        syms time
        % - g is the function that each of its elements multiplies one c_i 
        g = ((t - ti) ^ ord) * ((t - tf) ^ ord) * P;
        dgdt = diff(g, t);
        d2gdt2 = diff(g, t, 2);

        % - Construction Hermite interpolating polynomial -
        Ht = sym('Ht', [1 2 * ord]);
        Ht(1) = 1;
        for kk = 2:2*ord
            Ht(kk) = Ht(kk - 1) * t;
        end
        dnHdtn = sym('dnHdtn', [ord 2*ord]);
        for kk = 1:ord
            dnHdtn(kk, :) = diff(Ht, t, (kk - 1));
        end
        B = sym('B',[2*ord 1]);
        A = [subs(dnHdtn, t, ti)
            subs(dnHdtn, t, tf)];

        a = A \ B;
        H = Ht * a;
        dHdt = diff(H, t);
        d2Hdt2 = diff(dHdt, t);

        clear Ht kk dnHdtn A
        % - Time discretization and formation of the matrices
        dt = 1/Fs;                          % Sampling period (Original - sec)
        nt  = floor((tf - ti)/dt) + 1;      % Length of signal
        tt = (ti + (0:nt-1)*dt).';          % Time vector

        g0 = double(subs(g, tt));
        g1 = double(subs(dgdt, tt));
        g2 = double(subs(d2gdt2, tt));

        H0 = subs(H, tt);
        H1 = subs(dHdt, tt);
        H2 = subs(d2Hdt2, tt);

        MTRX_H = zeros([3, nt, 2*ord]);
        MTRX_H(1,:,:) = double(equationsToMatrix(H0, B));
        MTRX_H(2,:,:) = double(equationsToMatrix(H1, B));
        MTRX_H(3,:,:) = double(equationsToMatrix(H2, B));
    
        val = NaN([N/Nt, 1]);
        flag = NaN([N/Nt, 1]);
        if parallel == 1
            parfor_progress(N/Nt);
            parfor r = 1:N/Nt
                C0 = zeros(h, 1);
                options = optimset('Display', 'off');

                Herm = zeros([nt, 3*ndof]);
                for i = 1:ndof
                    for j = 1:3
                        Herm(:,3*(i-1)+j) = squeeze(MTRX_H(j,:,:)) * con(i,:,r).';
                    end
                end 
                % - Solve system of nonlinear eq. for calculating coeff. -
                [X, ~, flag(r)] = fsolve(@(C)nonlinsys(C, g0, g1, g2, Herm,...
                    dt, m0, c0, k0, e1, lambda), C0, options); %-- Change due to fractional

                if flag(r) > 0 
                    y = g0 * X + Herm(:,1);
                    y1 = g1 * X + Herm(:,2);
                    y2 = g2 * X + Herm(:,3);

                    lagrangian = (1/4).*pi.^(-1).*S0.^(-1).*(k0.*y+c0.*y1+m0.*y2+e1.*k0.*abs(y)).^2;

                    val(r) = trapz(tt, lagrangian);
                end  
                parfor_progress;
            end
            parfor_progress(0);
        else
            for r = 1:N/Nt
                r

                C0 = zeros(h, 1);
                options = optimset('Display', 'off');

                Herm = zeros([nt, 3*ndof]);
                for i = 1:ndof
                    for j = 1:3
                        Herm(:,3*(i-1)+j) = squeeze(MTRX_H(j,:,:)) * con(i,:,r).';
                    end
                end 
                % - Solve system of nonlinear eq. for calculating coeff. -
                [X, ~, flag(r)] = fsolve(@(C)nonlinsys(C, g0, g1, g2, Herm,...
                    dt, m0, c0, k0, e1, lambda), C0, options); %-- Change due to fractional

                if flag(r) > 0 
                    y = g0 * X + Herm(:,1);
                    y1 = g1 * X + Herm(:,2);
                    y2 = g2 * X + Herm(:,3);

                    lagrangian = (1/4).*pi.^(-1).*S0.^(-1).*(k0.*y+c0.*y1+m0.*y2+e1.*k0.*abs(y)).^2;

                    val(r) = trapz(tt, lagrangian);
                end            
            end
        end
        action = [action; val];
        exit_flag = [exit_flag; flag];
    end
    elapsed_time = toc;
end