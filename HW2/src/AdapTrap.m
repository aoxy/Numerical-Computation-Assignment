clear,clc
MS = 'MarkerSize'; ms = 20;

F = @(x) exp(-x.^2);
a = 0;
b = 1;
tol = 1e-6; % Preset tolerance
n = 1; % initial number of the subintervals
%%
h = (b-a) / n; % subinterval length
x = linspace(a, b, n+1); % points
% new value obtained by the composite trapezoidal rule
Tnew = h * ( F(a) / 2 + sum( F(x(2:end-1)) ) + F(b) / 2 ); 
%%
Told = Tnew + 1; % set an arbitrary value for the old value  
j = 1; 
while abs(Told-Tnew) > tol
    Told = Tnew; % keep the new quadrature value
    H = h * sum( F( a + (2*(1:n)-1)*h/2 ) );  % contribution from the midpoints
    Tnew = (Told + H) / 2;    % new value
    h = h / 2; % halve the subintervals
    n = 2*n;      % number of points in the new quadrature 
    plot(j, Tnew, '.', MS, ms) 
    hold on
    j = j+1; % number of the refinements
end

% print the results:
I = Tnew                                       
syms t
f = F(t);
Iexact = double(int(f, a, b))
AbsErr = abs(I - Iexact)
RelErr = AbsErr / abs(Iexact)
n