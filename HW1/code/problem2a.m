clear, clc, clf
LW = 'linewidth'; lw = 2;
nn = zeros(7, 1);
maxq = zeros(7, 1);

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
    %% 绘制图像
    figure(1);
    title('第一类边界条件样条插值效果');
    maxq(q - 5) = CubicSpline(x, F, h, M1, q); hold on
end

figure(2)
p3 = loglog(nn, maxq, 'b^-', 'MarkerFaceColor', 'b'); hold on
legend(p3, 'max error', 'location', 'se');
title('最大误差随n的log-log图');

function maxqq = CubicSpline(x, F, h, M, q)
LW = 'linewidth'; lw = 2;
n = size(x) - 1;
f = F(x);
areamax = zeros(2^q, 1);

for k = 1:n
    m = 4;
    xx = linspace(x(k), x(k + 1), m)';
    S = ((x(k + 1) - xx).^3 * M(k) + (xx - x(k)).^3 * M(k + 1)) / (6 * h(k)) + ...
        ((x(k + 1) - xx) * f(k) + (xx - x(k)) * f(k + 1)) / h(k) - ...
        h(k) * ((x(k + 1) - xx) * M(k) + (xx - x(k)) * M(k + 1)) / 6;
%     
%     p1 = plot(xx, F(xx), 'k', LW, lw); hold on
%     p2 = plot(xx, S, 'r', LW, lw); hold on
%     legend([p1, p2], 'exact', 'interpolant');
    error = abs(F(xx) - S);
    areamax(k) = max(error);
end

maxqq = max(areamax);
end
