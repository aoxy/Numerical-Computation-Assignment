clear, clc, clf
x = linspace(-1, 1, 6);
w = zeros(3, 1)';

for i = 1:3
    syms t;
    Fi = @(t) alpha_fun(t, x, i);
    % 使用课本130页下方的\alpha公式初始化w
    w(i) = int(Fi, t, -1, 1);
end

x = x(1:3);
dwdx = ones(6);

while (max(dwdx) > 1e-10)
    J = Jacobian(w, x);
    b = cal_b(w, x);
    dwdx = J \ b;
    dwdx = dwdx';
    w = w - dwdx(1:3);
    x = x - dwdx(4:6);
end
%% 输出w和x
w
x

function f = alpha_fun(t, x, i)
    n = length(x);
    fac1 = t - x;
    fac2 = x(i) - x;
    f = prod(fac1((1:n) ~= i))/prod(fac2((1:n) ~= i));
end
%% 第三题(b)中非线性方程组的Jacobian的表达式
function [J] = Jacobian(w,x)
    J = [1, 1, 1,0, 0, 0;
        x.^2,x .* w .* 2;
        x.^4,(x.^3) .* w .* 4;
        x.^6,(x.^5) .* w .* 6;
        x.^8,(x.^7) .* w .* 8;
        x.^10,(x.^9) .* w .* 10
        ];
end

function [b] = cal_b(w, x)
    b = zeros(6, 1);
    for i = 1:6
        b(i) = sum(w .* (x.^(2*i-2)))-1/(2*i-1);
    end
end
