clear, clc, clf
LW = 'linewidth'; lw = 2;

x = [-0.7 -0.5 0.25 0.75];
y = [0.99 1.21 2.57 4.23];
y1 = log(y);
%%
% y1 = ln(y) = ln(a) + bx
% aa(1) = ln(a), aa(2) = b
M = zeros (2, 2);
bb = zeros (2, 1);
xx = linspace (-1, 1, 1000);
%% 线性拟合
M(1, 1) = 4;

for i = 1:4
    M(1, 2) = M(1, 2) + x(i);
    M(2, 1) = M(2, 1) + x(i);
    M(2, 2) = M(2, 2) + (x(i))^2;
    bb(1) = bb(1) + y1(i);
    bb(2) = bb(2) + x(i) * y1(i);
end

aa = M \ bb;
a = exp(aa(1));
b = aa (2);
F = @(x) a * exp(b * x);
figure(1)
p1 = plot(x, y, 'o', LW, lw); hold on
plot(xx, F(xx), LW, lw);
h = legend('$$y_i$$', sprintf('$$y=%fe^{%fx}$$', a, b));
set(h, 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
%% 计算拟合函数的误差的2-范数
format long
error = abs(F(x) - y);
norm2 = sqrt(sum(error .* error))
%% 输出 norm2 = 0.006154650408178