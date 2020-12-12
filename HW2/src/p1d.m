clear, clc, clf
MS = 'MarkerSize'; ms = 8;
MC = 'MarkerFaceColor';

F = @(x) sin(cos(sin(cos(x)))); % 积分函数
a = -1; b = 1; % 积分区间[a,b]
tol = 1e-15; % 误差控制精度
n = 1; % 初始分点数

h = (b - a) / n; % 小区间长度
x = linspace(a, b, n + 1); % 取值点
% 新的复化梯形积分数值
Tnew = h * (F(a) / 2 + sum(F(x(2:end - 1))) + F(b) / 2);
T = zeros(100, 1); % 存储每个次迭代的复化梯形积分数值
Told = Tnew + 1; % 该Told无意义，用于进入下面的while循环
j = 1; % 迭代轮次
T(j) = Tnew; % 保存第1个轮的复化梯形积分数值

while abs(Told - Tnew) > tol
    Told = Tnew; % 保存上一轮迭代的复化梯形积分数值
    H = h * sum(F(a + (2 * (1:n) - 1) * h / 2));
    Tnew = (Told + H) / 2; % 计算当前轮的复化梯形积分数值
    T(j + 1) = Tnew;
    h = h / 2; % 小区间长度减半
    n = 2 * n; % 分点加密一倍
    j = j + 1; % 下一轮迭代
end

%% 输出计算结果
format long
I = Tnew%使用自动控制误差的复化梯形公式计算的积分值
syms t
f = F(t);
Iexact = double(int(f, a, b)); % 调库计算出的精确积分值
AbsErr = abs(I - Iexact)% 输出绝对误差
RelErr = AbsErr / abs(Iexact)% 输出相对误差
n% 输出最终分点数
semilogy((1:j), abs(T(1:j) - Iexact), 'ro-', MC, 'b', MS, ms);
xlabel('迭代次数');
ylabel('绝对误差');
