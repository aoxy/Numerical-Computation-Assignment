clear,clc, clf
MS = 'MarkerSize'; ms = 8;
MC = 'MarkerFaceColor';
%% 第1步
F = @(x) sin(cos(sin(cos(x)))); % 积分函数
a = -1; b = 1; % 积分区间[a,b]
e=1e-15; % 精度控制值
M=30; %最大循环次数
n=1;
h = b-a;
I = zeros(M); %存储每一轮迭代计算得到的积分值
%% 第2步
R= zeros(M,M);
R(1,1) = (F(a)+F(b))*h/2;
I(1) = R(1,1);
%% 第3步
for k=2 :M
    hk = h/2^(k-1);
    Fsum=sum(F(a + (2 * (1:2^(k-2)) - 1) * hk));
    R(k,1) = (R(k-1,1)+(2*hk)* Fsum)/2; % 2 * h_{k}=h_{k-1}
    for j=2:k
        R(k,j) = R(k,j-1)+(R(k,j-1)-R(k-1,j -1))/(4^(j-1) -1);
    end
    I(k) = R(k,k);
    n=n+1;
    if( abs(R(k,k )-R(k-1,k-1)) < e)
        Irichardson = R(k,k) % 第4步,输出最终的积分值
        break
    end
end
%% 输出计算结果
format long
n %迭代次数
semilogy((1:n), abs(I(1:n)-Irichardson), 'ro-', MC, 'b',MS,ms);
xlabel('迭代次数');
ylabel('误差');