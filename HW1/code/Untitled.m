clear, clc, clf
n=100;
lambda=((2^(n-1))/n) .*ones(n+1,1);
x=zeros(n+1,1);
for j=1:n+1
    x(j)=cos((j-1)*pi/n);
end
for j=2:n
    for k=1:n
        if j~=k
        lambda(j)=lambda(j)*(x(j)-x(k));
        end
    end
end
