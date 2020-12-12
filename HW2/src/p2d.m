clear, clc, clf
MC = 'MarkerFaceColor';

F = @(x) sin(3 * x.^2);
dF = @(x) 6 * x .* cos(3 * x.^2);
maxe = zeros(30, 1)';
maxc = zeros(30, 1)';
nvec = (1:2:59); % 依次令 n=1,3,5,...,57,59
j = 1;

for n = nvec
    %% 等距点
    xe = linspace(-1, 1, n + 1)';
    dfe = dF(xe); %精确导数值
    fe = F(xe); % 等距点上的函数值
    dpe = diffM(xe) * fe; %用微分矩阵计算的导数值
    maxe(j) = max(abs(dpe - dfe)); %误差绝对值最大值
    %% Chebyshev点
    xc = cos((0:n) * pi ./ n)'; % xc(j)中是x_{j-1}
    dfc = dF(xc);
    fc = F(xc);
    dpc = diffM(xc) * fc;
    maxc(j) = max(abs(dpc - dfc));
    j = j + 1;
end
%% 画图
figure(1)
plote = semilogy(nvec, maxe, 'b^-', MC, 'b'); hold on
plotc = semilogy(nvec, maxc, 'k*-', MC, 'k'); hold on
legend([plote, plotc], '等距点误差', 'Chebyshev点误差');
%% 计算微分矩阵D的函数
function D = diffM(x)
n = length(x);
D = zeros(n, n);
PI = zeros(n);
for k = 1:n
    PI(k) = prod(x(k) - x((1:n) ~= k));
end
for i = 1:n
    for j = 1:n
        if (i ~= j)
            D(i, j) = PI(i) / (PI(j) * (x(i) - x(j)));
        else
            D(i, j) = sum(1 ./ (x(j) - x((1:n) ~= j)));
        end
    end
end
end
