clear,clc,clf
A = [4,-2,4,2; -2,10,-2,-7; 4,-2,8,4; 2,-7,4,7];
b = [8; 2; 16; 6];
x = deEquations(A,b)
%% 运行的输出结果
% x =
%
%      1
%      2
%      1
%      2

%% 用Cholesky分解解方程组Ax=b
function [X]=deEquations(A,b)
N=length(A);
L=cholesky(A);
Lt=L';
X=zeros(N, 1);
Y=zeros(N, 1);
for j=1:N
    Y(j)=(b(j)-L(j,1:j-1)*Y(1:j-1))/L(j,j);
end
for k=N:-1:1
    X(k)=(Y(k)-Lt(k,k+1:N)*X(k+1:N))/Lt(k,k);
end
end

function [L]=cholesky(A)
N=length(A);
for i=1:N
    A(i,i)=sqrt(A(i,i)-A(i,1:i-1)*A(i,1:i-1)');
    if A(i,i)==0
        ME=MException('Zero Element', 'A(%d,%d) = 0',i,i);
        throw(ME);
    end
    for j=i+1:N
        A(j,i)=(A(j,i)-A(j,1:i-1)*A(i,1:i-1)')/A(i,i);
    end
end
L = tril(A); %取下三角部分
end
