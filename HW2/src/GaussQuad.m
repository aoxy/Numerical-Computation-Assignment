clear,clc
LW = 'LineWidth'; lw = 2;
MS = 'MarkerSize'; ms = 20;

F = @(x) x.^2.*cos(x);
N = 20;
nvec = 1:N; % number of the quadrature points
AbsErr = zeros(N, 1);
RelErr = zeros(N, 1);
syms x
Iexact = double(int(x^2*cos(x), x, -1, 1))
for n = nvec
    [x, w] = gauss(n);
    I = w*F(x);
    AbsErr(n) = abs(I - Iexact);
    RelErr(n) = AbsErr(n) / abs(Iexact);
end
semilogy(nvec, AbsErr, '.', MS, ms), hold on
semilogy(nvec, RelErr, 'o', LW, lw)
xlabel('number of quadrature points n')
ylabel('error')
legend('absolute Error', 'relative error')