% To solve y'=y-t^2+1 s.t. y(0) = 0.5 for 0<=t<=2 using forward Euler.
clear, clc, clf
LW = 'LineWidth'; lw = 2;
MS = 'MarkerSize'; ms = 18;

f = @(t, y) y-t^2+1;
h = 0.1;
T = 2;
N = T/h;

t = zeros(N+1, 1);
y = zeros(N+1, 1);
y(1) = 0.5;

for i = 1:N
    t(i+1) = i*h;
    y(i+1) = y(i) + h*f(t(i), y(i)); % forward Euler
end
plot(t, y, '.-', LW, lw, MS, ms), hold on

Y = @(t) (t+1).^2-.5*exp(t); % exact
tt = linspace(0, 2, 1000);
yy = Y(tt);
plot(tt, yy, 'k', LW, lw)

legend('Euler', 'exact', 'location', 'nw')
