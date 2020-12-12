% Comparison of finite difference schemes.
clear, clc, clf
MS = 'MarkerSize';
ms = 15;

f = @(x) sin(x); % the given function
fprime = @(x) cos(x); % the derivative
x = 1; % test point
H = 10.^-(1:0.5:3); % various size of h

err_fd = zeros(5, 1);
err_bd = zeros(5, 1);
err_cd = zeros(5, 1);

for k = 1:5
    h = H(k);
    
    D_fd = (f(x+h)-f(x))/h; % forward difference
    D_bd = (f(x)-f(x-h))/h; % forward difference
    D_cd = (f(x+h)-f(x-h))/(2*h); % center difference
    
    err_fd(k) = abs(D_fd-fprime(x));
    err_bd(k) = abs(D_bd-fprime(x));
    err_cd(k) = abs(D_cd-fprime(x));    
end

loglog(H.', err_fd, '.-k', MS, ms), hold on
loglog(H.', err_bd, '.-r', MS, ms)
loglog(H.', err_cd, '.-b', MS, ms)

legend('forward difference (first order)', ...
    'backward difference (first order)', ...
    'central difference (second order)', 'location', 'se')

xlabel('step size h')
ylabel('error at point x')
title('Comparison of finite difference schemes')