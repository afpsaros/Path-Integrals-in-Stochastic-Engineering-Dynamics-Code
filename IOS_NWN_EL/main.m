clc
clear all
close all

% System parameters

system_params.m0 = 1;
system_params.c0 = 0.1;
system_params.k0 = 1;
system_params.S0 = 0.0637;
system_params.e1 = 0;
system_params.e2 = 0;
system_params.p0 = 1;
system_params.q0 = 0.1;
system_params.r0 = 1;
%%
points = [5 5 5 5 1];
bounds = [-2 2 -2 2 -2 2 -2 2 1 1];

x1 = linspace(bounds(1), bounds(2), points(1));
x2 = linspace(bounds(3), bounds(4), points(2));
x3 = linspace(bounds(5), bounds(6), points(3));
x4 = linspace(bounds(7), bounds(8), points(4));
times = linspace(bounds(9), bounds(10), points(5));

[X1, X2, X3, X4, TIMES] = ndgrid(x1, x2, x3, x4, times); 

eval_params.grid = [X1(:) X2(:) X3(:) X4(:) TIMES(:)];
eval_params.N = prod(points);
eval_params.Nt = points(5);
eval_params.ti = 0;
%%
algo_params.numt = 20;
algo_params.numinit = 20;
algo_params.maxiter = 10000;
algo_params.tol = 1e-4;
algo_params.parallel = 0;
algo_params.setopt = 1;
%%
[action, elapsed_time] = ELminim(system_params, eval_params, algo_params);    

%%
pdft = reshape(exp(-action), points(1), points(2), points(3), points(4), points(5));

for i = 1:points(5)
    pdf = reshape(pdft(:, :, :, :, i), points(1), points(2), points(3), points(4));
    f123 = trapz(x4, pdf, 4);
    f12 = trapz(x3, f123, 3);
    f1 = trapz(x2, f12, 2);
    f1 = f1 / trapz(x1, f1);
    f1t(:, i) = f1.';
end

for i = 1:points(5)
    pdf = reshape(pdft(:, :, :, :, i), points(1), points(2), points(3), points(4));
    f123 = trapz(x4, pdf, 4);
    f12 = trapz(x3, f123, 3);
    f2 = trapz(x1, f12, 1);
    f2 = f2 / trapz(x2, f2);
    f2t(:, i) = f2.';
end
%%
figure;
plot(x1, f1t)
xlabel('x')
ylabel('PDF')

figure;
plot(x2, f2t)
xlabel('v')
ylabel('PDF')
