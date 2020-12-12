clear, clc, clf
MS = 'MarkerSize'; ms = 8;
MC = 'MarkerFaceColor';

F = @(x) sin(cos(sin(cos(x)))); % 积分函数
a = -1; b = 1; % 积分区间[a,b]
e = 1e-15; % 精度控制值
M = 30; %最大循环次数
I = zeros(M);
[x, w] = gauss(1);
I(1) = w * F(x);

for n = 2:M
    [x, w] = gauss(n);
    I(n) = w * F(x);
    
    if (abs(I(n) - I(n - 1)) < e)
        Igauss = I(n)% 输出最终的积分值
        break
    end
    
end


n%迭代次数
semilogy((1:n), abs(I(1:n) - Igauss), 'ro-', MC, 'b', MS, ms);
xlabel('迭代次数');
ylabel('误差');

% GAUSS  nodes x (Legendre points) and weights w
%        for Gauss quadrature
function [x, w] = gauss(N)
beta = .5 ./ sqrt(1 - (2 * (1:N - 1)).^(-2));
T = diag(beta, 1) + diag(beta, -1);
[V, D] = eig(T);
x = diag(D); [x, i] = sort(x);
w = 2 * V(1, i).^2;
end
