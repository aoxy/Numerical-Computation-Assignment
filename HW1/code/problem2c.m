clear, clc, clf
LW = 'linewidth'; lw = 2;

nn = zeros(7, 1);
maxq1 = zeros(7, 1);
maxq2 = zeros(7, 1);
maxq3 = zeros(7, 1);

for q = 6:12
    nn(q - 5) = 2^q;
    n = 2^q;
    x = linspace(-1, 1, n + 1)';
    F = @(x)exp (3 .* cos(pi .* x));
    f = F(x);
    h = diff(x);
    df = diff(f);
    lambda = h(2:n) ./ (h(2:n) + h(1:n - 1));
    d = 6 * (df(2:n) ./ h(2:n) - df(1:n - 1) ./ h(1:n - 1)) ./ (h(2:n) + h(1:n - 1));
    mu = 1 - lambda;
    
    %% 第一类边界条件
    M0 = 0;
    Mn = 0;
    A1 = diag(2 * ones(n - 1, 1)) + diag(lambda(1:n - 2), 1) + diag(mu(2:n - 1), -1);
    D1 = [d(1) - mu(1) * M0; d(2:n - 2); d(n - 1) - lambda(n - 1) * Mn];
    M1 = A1 \ D1;
    M1 = [M0; M1; Mn];
    %% 第二类边界条件
    m0 = 0;
    mn = 0;
    lambda2 = [1; lambda];
    mu2 = [mu; 1];
    d0 = 6 * (df(1) / h(1) - m0) / h(1);
    dn = 6 * (mn - df(n) / h(n)) / h(n);
    D2 = [d0; d; dn];
    A2 = diag(2 * ones(n + 1, 1)) + diag(lambda2, 1) + diag(mu2, -1);
    M2 = A2 \ D2;
    %% 第三类边界条件
    lambda0 = h(1) / (h(1) + h(n));
    lambda3 = [lambda0; lambda(1:n - 2)];
    mu0 = 1 - lambda0;
    d0 = 6 * (df(1) ./ h(1) - df(n) ./ h(n)) / (h(1) + h(n));
    D3 = [d0; d];
    A3 = diag(2 * ones(n, 1)) + diag(lambda3, 1) + diag(mu, -1);
    A3(1, n) = mu0;
    A3(n, 1) = lambda(n - 1);
    M3 = A3 \ D3;
    M3 = [M3; M3(1)];
    %% 求最大误差
    maxq1(q - 5) = CubicSpline(x, F, h, M1, q);
    maxq2(q - 5) = CubicSpline(x, F, h, M2, q);
    maxq3(q - 5) = CubicSpline(x, F, h, M3, q);
end

%% 绘制图像
figure(1)
p3 = loglog(nn, maxq1, 'b^-', 'MarkerFaceColor', 'b'); hold on
p4 = loglog(nn, maxq2, 'ro-'); hold on
p5 = loglog(nn, maxq3, 'k*-', 'MarkerFaceColor', 'k'); hold on
legend([p3, p4, p5], '第一类边界条件', '第二类边界条件', '第三类边界条件')
title('三种边界条件下最大误差随n的log-log图')

function maxqq = CubicSpline(x, F, h, M, q)
n = size(x) - 1;
f = F(x);
areamax = zeros(2^q, 1);

for k = 1:n
    m = 4;
    xx = linspace(x(k), x(k + 1), m)';
    S = ((x(k + 1) - xx).^3 * M(k) + (xx - x(k)).^3 * M(k + 1)) / (6 * h(k)) + ...
        ((x(k + 1) - xx) * f(k) + (xx - x(k)) * f(k + 1)) / h(k) - ...
        h(k) * ((x(k + 1) - xx) * M(k) + (xx - x(k)) * M(k + 1)) / 6;
    error = abs(F(xx) - S);
    areamax(k) = max(error);
end

maxqq = max(areamax);
end
