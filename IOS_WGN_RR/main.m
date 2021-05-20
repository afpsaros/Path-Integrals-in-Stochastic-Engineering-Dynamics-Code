clc
clear all
close all

% System parameters

system_params.m0 = 1;
system_params.c0 = 0.1;
system_params.k0 = 1;
system_params.S0 = 0.0637;
system_params.e1 = 0.1;
system_params.e2 = 0.1;
system_params.lambda = 0;
system_params.ord = 2;
system_params.ndof = 1;
%%
points = [5 5 2];
bounds = [-2 2 -2 2 1 2];

x1 = linspace(bounds(1), bounds(2), points(1));
x2 = linspace(bounds(3), bounds(4), points(2));
times = linspace(bounds(5), bounds(6), points(3));

[X1, X2, TIMES] = ndgrid(x1, x2, times); 

eval_params.grid = [X1(:) X2(:) TIMES(:)];
eval_params.N = prod(points);
eval_params.Nt = points(3);
eval_params.ti = 0;
%%
algo_params.h = 5;
algo_params.Fs = 100;
algo_params.parallel = 1;
%%
[action, exit_flag, elapsed_time] = RRminim(system_params, eval_params, algo_params);
pdft = reshape(exp(-action), points(1), points(2), points(3));

for i = 1:points(3)
    pdf = reshape(pdft(:, :, i), points(1), points(2));
    f1 = trapz(x2, pdf, 2);
    f1 = f1 / trapz(x1, f1);
    f1t(:, i) = f1.';
end

for i = 1:points(3)
    pdf = reshape(pdft(:, :, i), points(1), points(2));
    f2 = trapz(x1, pdf, 1);
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
