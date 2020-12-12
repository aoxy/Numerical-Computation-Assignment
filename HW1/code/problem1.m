clear, clc, clf
LW = 'linewidth'; lw = 2;

n = 5000;
x = zeros(n + 1, 1);
m = 10000;
xx = linspace(-1, 1, m)';
F = @(x)tanh(20 * sin(12 .* x)) + 0.02 * exp(3 .* x) .* sin(300 .* x);
f = F(x);
p1 = zeros(m, 1);
p2 = zeros(m, 1);
p = zeros(m, 1);
%% 基于Chebyshev点的第二形式的重心插值公式
for j = 1:n + 1
    x(j) = cos((j - 1) * pi / n);
end

p1 = 0.5 * (F(1) ./ (xx -1) + ((-1)^n) * F(-1) ./ (xx + 1));
p2 = 0.5 * (1 ./ (xx -1) + ((-1)^n) ./ (xx + 1));

for k = 2:n
    p1 = p1 + ((-1)^(k - 1)) * F(x(k)) ./ (xx - x(k));
    p2 = p2 + ((-1)^(k - 1)) ./ (xx - x(k));
end

p = p1 ./ p2;
%% 被插值函数图像
figure(1)
plot(xx, F(xx), 'r', LW, 4), hold on
legend('exact', 'location', 'nw')
%% 误差图像
figure(2)
plot(2)
semilogy(xx, abs(F(xx) - p), 'k', LW, lw), hold on
legend('error', 'location', 'se')