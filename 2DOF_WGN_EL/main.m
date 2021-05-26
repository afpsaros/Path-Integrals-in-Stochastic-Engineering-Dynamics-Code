clc
clear all
close all

% Number of data points
N = 100;

dim = 5;
hal = haltonset(dim);

bounds_x1 = [-4 4];
bounds_x2 = [-4 4];
bounds_x3 = [-4 4];
bounds_x4 = [-4 4];
bounds_t = [1 1];
bounds_prod = repmat([bounds_x1(2)-bounds_x1(1) bounds_x2(2)-bounds_x2(1) bounds_x3(2)-bounds_x3(1) bounds_x4(2)-bounds_x4(1) bounds_t(2)-bounds_t(1)], N, 1);
bounds_shift = repmat([bounds_x1(1) bounds_x2(1) bounds_x3(1) bounds_x4(1) bounds_t(1)], N, 1);

dsites = net(hal, N) .* bounds_prod + bounds_shift;
ctrs = dsites;
%
m = 1;
c = 0.35;
k = 0.5;
e1 = 0.1;
e2 = 0.1;

S0 = 0.1;

t = dsites(:, 5).';
grid = dsites(:, 1:4);
%%
neval = 5; nevalt = 1; M = neval^4 * nevalt; 
x1 = linspace(bounds_x1(1), bounds_x1(2), neval);
x2 = linspace(bounds_x2(1), bounds_x2(2), neval);
x3 = linspace(bounds_x3(1), bounds_x3(2), neval);
x4 = linspace(bounds_x4(1), bounds_x4(2), neval);
t = linspace(bounds_t(1), bounds_t(2), nevalt);
[X1, X2, X3, X4, T] = ndgrid(x1, x2, x3, x4, t); epoints = [X1(:) X2(:) X3(:) X4(:) T(:)];
N = M;
twpi = epoints(:, 5).';
grid = epoints(:, 1:4);
%%
numt = 100;
numinit = 20;
maxiter = 50;
tol = 10^-3;
parallel = 1;
setopt = 1;

tic
val = ELminim_tdof(m, c, k, e1, e2, twpi, N, grid, S0, numt, numinit, ...
        maxiter, tol, setopt, parallel);
toc
%%
val = inpaintn(val,200); %%%%%%%%%%%%%%%%%%%%%%%%
%%
x1f = linspace(bounds_x1(1), bounds_x1(2), 101);
x2f = linspace(bounds_x2(1), bounds_x2(2), 101);
Pf = reshape(exp(-val), neval, neval, neval, neval, nevalt);
for i = 1:nevalt
    pdf = reshape(Pf(:, :, :, :, i), neval, neval, neval, neval);
    f123 = trapz(x4, pdf, 4);
    f12 = trapz(x3, f123, 3);
    f1 = trapz(x2, f12, 2);
    f1 = f1 / trapz(x1, f1);
    f1 = interp1(x1, f1, x1f, 'linear');
    f1t(:, i)=f1.';
end

for i = 1:nevalt
    pdf = reshape(Pf(:, :, :, :, i), neval, neval, neval, neval);
    f234 = reshape(trapz(x1, pdf, 1), neval, neval, neval);
    f23 = trapz(x4, f234, 3);
    f2 = trapz(x3, f23, 2);
    f2 = f2 / trapz(x2, f2);
    f2 = interp1(x2, f2, x2f, 'linear');
    f2t(:, i)=f2.';
end

for i = 1:nevalt
    pdf = reshape(Pf(:, :, :, :, i), neval, neval, neval, neval);
    f234 = reshape(trapz(x1, pdf, 1), neval, neval, neval);
    f34 = reshape(trapz(x2, f234, 1), neval, neval);
    f3 = trapz(x4, f34, 2);
    f3 = f3 / trapz(x3, f3);
    f3 = interp1(x3, f3, x1f, 'linear');
    f3t(:, i)=f3.';
end

for i = 1:nevalt
    pdf = reshape(Pf(:, :, :, :, i), neval, neval, neval, neval);
    f234 = reshape(trapz(x1, pdf, 1), neval, neval, neval);
    f34 = reshape(trapz(x2, f234, 1), neval, neval);
    f4 = reshape(trapz(x3, f34, 1), neval,1);
    f4 = f4 / trapz(x4, f4);
    f4 = interp1(x4, f4, x2f, 'linear');
    f4t(:, i)=f4.';
end
%%
figure;
plot(x1f, f3)
figure;
plot(x2f, f4)