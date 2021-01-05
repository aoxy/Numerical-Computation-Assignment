n=2;
A=eye(n);
A(1,1)=1;
A(1,2)=2;
A(2,1)=3;
A(2,2)=4;
B =[1, 2, 3;2, 3, 4;3, 4, 5];
C=[1, -2, -2, -3;3, -9, 0, -9;-1, 2, 4, 7;-3, -6, 26, 2];
dDool(C)
%LU分解法求解Ax=b，假定A矩阵可进行LU分解以及对角线元素均不为0
function dDool(A)
n=length(A);A(2:n,1)=A(2:n,1)/A(1,1);
for t=2:n-1  %进行LU分解
    A(t,t:n)=A(t,t:n)-A(t,1:t-1)*A(1:t-1,t:n);
    t
    t:n
    A(t+1:n,t)=(A(t+1:n,t)-A(t+1:n,1:t-1)*A(1:t-1,t))/A(t,t);
end
A(n,n)=A(n,n)-A(n,1:n-1)*A(1:n-1,n);
L=tril(A,-1)+eye(n);U=triu(A); %矩阵A分解出L和U
L
U
A

end