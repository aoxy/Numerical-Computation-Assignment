clear, clc, clf
LW = 'linewidth'; lw = 2;

n = 15;
x = linspace(-1, 1, n + 1)';
m = 1001;
xx = linspace(-1, 1, m)';
F = @sin;
f = F(x);
P = @cos;
p = P(xx);
pp = zeros(m, 1);

for k = 1:n
    l = ones(m, 1);
    sum = zeros(m, 1); % 公式(1)中大括号内的求和部分
    
    for j = 1:k - 1
        l = l .* (xx - x(j)) / (x(k) - x(j));
        sum = sum +1 ./ (xx - x(j));
    end
    
    for j = k + 1:n
        l = l .* (xx - x(j)) / (x(k) - x(j));
        sum = sum +1 ./ (xx - x(j));
    end
    
    pp = pp + f(k) * l .* sum;
end

R = abs(pp - p); %与真实解的误差绝对值
figure (1)
plot(xx, P(xx), 'k', LW, lw), hold on
plot(xx, pp, 'b', LW, lw);
xlabel('x');
ylabel('p(x)或cos(x)');
legend('exact', 'interpolant', 'location', 'nw');
figure (2)
semilogy(xx, R, 'k', LW, lw);
xlabel('x');
ylabel('error');
