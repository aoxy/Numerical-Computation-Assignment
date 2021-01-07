clear,clc,clf

[w,x]=gauss_w_x(6)


function [w,x]=gauss_w_x(n)
halfn=ceil(n/2);
w=zeros(halfn,1)';
x = linspace(-1, 1, n);
syms t;
for i=1:halfn
    Fi=@(t) alpha_fun(t,x,i);
    w(i)=int(Fi,t,-1,1);
end
x=x(1:floor(n/2));
dx=ones(halfn);
dw=ones(halfn);
while(max([dx,dw]) > 1e-5)
    J = Jacobian(x,w,n);
    b = f(x,w,n);
    [dx,dw] =grad(J,b);%J\b;
    x = x + dx;
    w = w + dw;
end
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

function [J] = Jacobian(x ,w,n)
J=zeros(n,n);
m=length(w);
if mod(n,2)
    factor=ones(n,1)';
    factor(m)=1/2;
    xp=[x,0];
    for i=1:n
        J(i,:) = [(xp .^ (2*i-3)) .* w .* (2*i-2) , x .^ (2*i-2)].*factor;
    end
    J
else
    for i=1:n
        J(i,:) = [(x .^ (2*i-3)) .* w .* (2*i-2) , x .^ (2*i-2)];
    end
end

%     x_p=x(1:n-1);
%     w_p=w(1:n-1);
%     %J(1,:) = [0,0,0,1,1,1];
%     for i=1:n
%         J(1,:) = [(x_p .^ (2*i-3)) .* w_p .* (2*i-2) ,, x_p .^ (2*i-2),];
%     end

%     J(1,1:6) = [0,0,0,1,1,1];
%     J(2,1:6) = [x .* w .* 2 , x .^ 2];
%     J(3,1:6) = [(x .^ 3) .* w .* 4 , x .^ 4];
%     J(4,1:6) = [(x .^ 5) .* w .* 6 , x .^ 6];
%     J(5,1:6) = [(x .^ 7) .* w .* 8 , x .^ 8];
%     J(6,1:6) = [(x .^ 9) .* w .* 10 , x .^ 10];
end

function [b] = f(x,w,n)
b = zeros(n,1);
if mod(n,2)
    xp=[x,0];
    for i=1:n
        b(i)=sum(w .* (xp .^ (2*i-2))) - 1/(2*i-1);
    end
else
for i=1:n
    b(i)=sum(w .* (x .^ (2*i-2))) - 1/(2*i-1);
end
end
% b(1) = sum(w) - 1;
% b(2) = sum(w .* (x .^ 2)) - 1/3;
% b(3) = sum(w .* (x .^ 4)) - 1/5;
% b(4) = sum(w .* (x .^ 6)) - 1/7;
% b(5) = sum(w .* (x .^ 8)) - 1/9;
% b(6) = sum(w .* (x .^ 10)) - 1/11;

b = b .* (-1);  % b直接用作方程组求解
end