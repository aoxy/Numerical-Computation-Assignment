clear,clc
    MS = 'MarkerSize'; ms = 20;
    
    F = @(x) sin(cos(sin(cos(x))));
    a = -1;
    b = 1;
    T=zeros(20,1) ; % 每轮计算的积分值
    n = 1; % initial number of the subintervals
    %%
    h = (b-a) / n; % subinterval length
    x = linspace(a, b, n+1); % points
    % new value obtained by the composite trapezoidal rule
    Tnew = h * ( F(a) / 2 + sum( F(x(2:end-1)) ) + F(b) / 2 ); 
    %%
    Told = Tnew + 1; % set an arbitrary value for the old value  
    j = 1; 
    N=25;
    while j<N      %算20轮误差，看看到后面误差是否会趋于稳定
        Told = Tnew; % keep the new quadrature value
        H = h * sum( F( a + (2*(1:n)-1)*h/2 ) );  % contribution from the midpoints
        Tnew = (Told + H) / 2;    % new value
        h = h / 2; % halve the subintervals
        n = 2*n;      % number of points in the new quadrature
        T(j)=Tnew; 
        figure(1)
        semilogy(j, abs(Told-Tnew)/3, '.', MS, ms)
        title('每轮计算出的递归误差：') 
        hold on
        j = j+1; % number of the refinements
    end
    %%得到精确值之后画出误差图：
    last=T(N-1)*ones(N-1,1);
    x_j=linspace(1,N-1,N-1);
    figure(2)
    semilogy(x_j,abs(T-last),'.', MS, ms)
    title('最终的log-log图：')
    xlabel('轮数')
    ylabel('每轮的误差')
    hold on

    %输出结果：
    I=Tnew
    syms t
    f=F(t)
    Iexact=double(int(f,a,b))
