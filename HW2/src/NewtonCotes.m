clear, clc

F = @(x) 2*sin(x) + exp(x.^2-1);
a = 0;
b = 1;
n = 4; % the degree of the underlying interpolant
%%
x = linspace(a, b, n+1)';
syms t
for i = 0:n
    L = 1;
    for j = 0:n
        if j ~= i
            L = L * (t-j)/(i-j); % construction of the Lagrange polynomial
        end
    end
    % symbolic calculation of the quadrature weights divided by the length 
    % of the interval:
    c(i+1) = int(L, 0, n) / n;
end
c % print
c = double(c); % conversion to floating point numbers
I = (b-a) * c * F(x)       % print the result
Iexact = double(int(F(t), a, b)) % exact value
AbsErr = abs(I - Iexact) % absolute error
RelErr = AbsErr / abs(Iexact) % relative error