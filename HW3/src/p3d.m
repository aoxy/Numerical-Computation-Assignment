clear,clc,clf
%% 没写完的垃圾代码2021年1月7日
[w,x]=gauss_w_x(5)

function [w,x]=gauss_w_x(n)
halfn=ceil(n/2);
w=zeros(n,1)';
x=che(1:n,n);
syms t;
for i=1:halfn
    Fi=@(t) alpha_fun(t,x,i);
    w(i)=int(Fi,t,-1,1);
end
dx=ones(halfn);
dw=ones(halfn);
while(max([dx,dw]) > 1e-5)

    J = Jacobian(x,w,n);
    b = f(x,w,n);
        mid=fix((n+1)/2);
    x(mid+1:n)=x(mid+1:n)+dxdw(1:n-mid);
    w(n-mid:n)=w(n-mid:n)+dwdw(n-mid:n);
for i=1 : n/2
x(i)=-x(n-i+1);
w(i)=w(n-i+1);
end
    dxdw=J^(-1)*b;
    %[dx,dw] =grad(J,b);
    %x = x - dx;
    %w = w - dw;
end
end
%% Chebyshev点
function [xj]=che(j,n)
    xj = cos(j*pi/n);
end
function [dx,dw] =grad(J,b)
n=size(J,1);
fn=floor(n/2);
dxdw=J\b;
dxdw=dxdw';
dx=dxdw(1:fn);
dw=dxdw(fn+1:n);
end
function f=alpha_fun(t,x,i)
n=length(x);
fac1=t-x;
fac2=x(i)-x;
f=prod(fac1((1:n) ~= i))/prod(fac2((1:n) ~= i));
end
function [J]= Jacobian(x,w,n)
J=zeros(n,n);
for i=1 :n
for j=1 : n
k=2*(i-1);
mid=fix(n/2);
l=fix((n+1)/2);
if j<=mid
if k==0
J(i,j)=0;
else
J(i,j)=k*w(j+l)*(x(j+l)^(k-1));
end
else
J(i,j)=x(j)^k;
end
end
end
% k=mod(n,2)
if mod(n,2)
J(1,(n+1)/2)=1/2;
end
end
function [J] = Jacobian1(x ,w,n)
J=zeros(n,n);
m=length(w);
if mod(n,2)
    factor=ones(n,1)';
    factor(m)=1/2;
    xp=[x,0];
    m1zeros=zeros(m,1)';
    m2ones=ones(length(x),1)';
    J(1,:) = [m1zeros , m2ones].*factor;
    for i=2:n
        J(i,:) = [(xp .^ (2*i-3)) .* w .* (2*i-2) , x .^ (2*i-2)].*factor;
    end
    J
else
    for i=1:n
        J(i,:) = [(x .^ (2*i-3)) .* w .* (2*i-2) , x .^ (2*i-2)];
    end
end
end
function [F]= f(x,w,n)
xf=zeros(n,1);
for i=1 : n
k=2*(i-1);
xf(i)=1/(k+1);
for j=1: (1+n)/2
xf(i)=xf(i)-w(j)*(-x(j))^k;
end
end
F=xf;
if mod(n,2)
F(1)=0;
end
end
function [b] = f1(x,w,n)
b = zeros(n,1);
if mod(n,2)
    w=w(1:end-1);
    b(1)=sum(w) - 1;
    for i=2:n
        b(i)=sum(w .* (x .^ (2*i-2))) - 1/(2*i-1);
    end
else
    for i=1:n
        b(i)=sum(w .* (x .^ (2*i-2))) - 1/(2*i-1);
    end
end
end