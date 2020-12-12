clear, clc, clf
LW = 'linewidth'; lw = 2;
a = 0; b = 2; % 求解区间[0, 2]
y = @(x) (x.^2) .* exp(-5 .* x) / 2; % 推导出的精确解
nvec = 101:10001; % 取样点数从101到1001
hvec = (b - a) ./ (nvec - 1);
maxe = zeros(length(nvec), 1)';
iter = 1; % 计数
for n = nvec
    h = (b - a) / (n - 1); %步长
    x = linspace(a, b, n)'; % 取样点
    g = x .* exp(-5 .* x);
    yn = zeros(n, 1); % 存储计算出的函数在取样点处的数值解
    yn(1) = 0; %初值条件
    f = @(x, y) x * exp(-5 * x) - 5 * y;
    %% 使用三阶Runge-Kutta方法起步，计算前两项即可
    k1 = zeros(2, 1);
    k2 = zeros(2, 1);
    k3 = zeros(2, 1);
    for i = 1:2
        k1(i) = f(x(i), yn(i));
        k2(i) = f(x(i) + h / 2, yn(i) + (h / 2) * k1(i));
        k3(i) = f(x(i) + h, yn(i) - h * k1(i) + 2 * h * k2(i));
        yn(i + 1) = yn(i) + (h / 6) * (k1(i) + 4 * k2(i) + k3(i));
    end
    %% 使用4阶的线性多步法从第3项开始计算
    for i = 3:n - 1
        yn(i + 1) = (yn(i - 1) + (h / 3) * (g(i - 1) - 5 * yn(i - 1) + 4 * g(i) - 20 * yn(i) + g(i + 1))) / (1 + 5 * h / 3);
    end
    ya = y(x); % 精确解在取样点处的取值
    maxe(iter) = max(abs(ya - yn));
    iter = iter + 1;
end

figure(1)
pn = loglog(hvec, maxe, 'k', LW, lw); hold on
xlabel('h')
ylabel('max error')
set(gca, 'yaxislocation', 'right');
