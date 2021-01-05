clear, clc, clf
%% 构造A和b
D = orth(rand(5, 5));
B = diag([1, 2, 3, 4, 5]);
A = D \ B * D;
b = rand(5, 1);
xexact = A \ b; % x精确解

%% 使用多个omega迭代求解
lambda1 = 1;
lambdan = 5;
omegab=2/(lambda1+lambdan);
eponch=1001;
omegalist = linspace(0.001, 0.399, eponch); 
count=zeros(eponch,1)';
for i = 1:eponch
    count(i)=Richardson(A,b,xexact,omegalist(i));
end

%% 输出结果
[mink,idx]=min(count);
omegalist(idx);
fp=fopen('p2c_out.txt','w');
fprintf (fp,'bestomega = %f\tk=%d\n',omegalist(idx),mink);
fprintf (fp,'omegab = %f\tk=%d\n',omegab,Richardson(A,b,xexact,omegab));
semilogy(omegalist, count, 'k.-');
xlabel('\omega');
ylabel('迭代次数');

%% 返回Richardson迭代方法的迭代次数
function [k]=Richardson(A,b,xexact,omega)
    G = eye(5) - omega * A;
    x = zeros(5, 1);
    k = 0;
    while norm(x - xexact) > 1e-13
        x = G * x + omega * b;
        k = k + 1;
    end
end