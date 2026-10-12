%% Exercise 3 SDE
% Fixed-step explicit-explicit SDE solver
% Fixed-step implicit-explicit SDE solver
% Test both on Geometric Brownian Motion
% Estimate numerical order from strong error at final time

clear
clc
close all

%% Geometric Brownian Motion setup from Lecture 03A
lambda = 0.15;
sigma = 0.15;
x0 = 1.0;

p = [lambda; sigma];

Tfinal = 10;
N = 1000;
Ns = 2000;
seed = 100;

%% Generate Wiener realizations
[W,T,dW] = ScalarStdWienerProcess(Tfinal,N,Ns,seed);

%% Solve with explicit-explicit Euler-Maruyama
Xexplicit = zeros(Ns,N+1);

for i = 1:Ns
    Xexplicit(i,:) = SDEeulerExplicitExplicit( ...
        @GeometricBrownianDrift, ...
        @GeometricBrownianDiffusion, ...
        T, ...
        x0, ...
        W(i,:), ...
        p);
end

%% Solve with implicit-explicit Euler
Ximplicit = zeros(Ns,N+1);

for i = 1:Ns
    Ximplicit(i,:) = SDEeulerImplicitExplicit( ...
        @GeometricBrownianDriftWithJacobian, ...
        @GeometricBrownianDiffusion, ...
        T, ...
        x0, ...
        dW(i,:), ...
        p);
end

%% Exact Geometric Brownian Motion solution on the same Wiener paths
Xexact = zeros(Ns,N+1);

for i = 1:Ns
    Xexact(i,:) = x0*exp( ...
        (lambda - 0.5*sigma^2)*T ...
        + sigma*W(i,:));
end

%% Plot a few realizations
nPlot = 8;

figure

subplot(3,1,1)
plot(T,Xexact(1:nPlot,:)','LineWidth',1.0)
xlabel('t')
ylabel('x(t)')
title('Exact Geometric Brownian Motion')
grid on

subplot(3,1,2)
plot(T,Xexplicit(1:nPlot,:)','LineWidth',1.0)
xlabel('t')
ylabel('x(t)')
title('Explicit-Explicit Euler-Maruyama')
grid on

subplot(3,1,3)
plot(T,Ximplicit(1:nPlot,:)','LineWidth',1.0)
xlabel('t')
ylabel('x(t)')
title('Implicit-Explicit Euler')
grid on

%% Plot one realization against exact solution
figure

plot(T,Xexact(1,:),'k','LineWidth',1.5)
hold on
plot(T,Xexplicit(1,:),'--','LineWidth',1.2)
plot(T,Ximplicit(1,:),':','LineWidth',1.2)

xlabel('t')
ylabel('x(t)')
legend('Exact','Explicit-explicit','Implicit-explicit','Location','best')
title('Same Wiener realization')
grid on

%% Final-state distributions
figure

histogram(Xexact(:,end),40,'Normalization','pdf')
hold on
histogram(Xexplicit(:,end),40,'Normalization','pdf')
histogram(Ximplicit(:,end),40,'Normalization','pdf')

xlabel('x(T)')
ylabel('Probability density')
legend('Exact','Explicit-explicit','Implicit-explicit','Location','best')
title('Distribution of final state')
grid on

%% Strong error at final time for the current step size
strongErrExplicit = mean(abs(Xexact(:,end)-Xexplicit(:,end)));
strongErrImplicit = mean(abs(Xexact(:,end)-Ximplicit(:,end)));

fprintf('Single-grid strong errors at T = %.2f\n',Tfinal)
fprintf('Explicit-explicit: %.6e\n',strongErrExplicit)
fprintf('Implicit-explicit: %.6e\n',strongErrImplicit)

%% Determine numerical order
% Lecture 03A defines order p by
% E{|x(T)-x_N|} = O(dt^p)
%
% We use the exact GBM solution and the SAME underlying Wiener paths
% for all step sizes, then fit
%
% log(error) = p*log(dt) + constant

Nfine = 6400;
NsOrder = 2000;
seedOrder = 321;

[Wfine,Tfine,dWfine] = ScalarStdWienerProcess( ...
    Tfinal,Nfine,NsOrder,seedOrder);

Nlist = [100 200 400 800 1600 3200];
dtList = Tfinal./Nlist;

errExplicit = zeros(size(Nlist));
errImplicit = zeros(size(Nlist));

for j = 1:length(Nlist)

    Nj = Nlist(j);
    block = Nfine/Nj;

    if abs(block-round(block)) > 0
        error('Nfine must be divisible by every element of Nlist.')
    end

    block = round(block);

    dWcoarse = zeros(NsOrder,Nj);

    for k = 1:Nj
        idx = (k-1)*block + (1:block);
        dWcoarse(:,k) = sum(dWfine(:,idx),2);
    end

    Wcoarse = [zeros(NsOrder,1) cumsum(dWcoarse,2)];
    Tcoarse = linspace(0,Tfinal,Nj+1);

    XexpFinal = zeros(NsOrder,1);
    XimpFinal = zeros(NsOrder,1);
    XexactFinal = zeros(NsOrder,1);

    for i = 1:NsOrder

        Xexp = SDEeulerExplicitExplicit( ...
            @GeometricBrownianDrift, ...
            @GeometricBrownianDiffusion, ...
            Tcoarse, ...
            x0, ...
            Wcoarse(i,:), ...
            p);

        Ximp = SDEeulerImplicitExplicit( ...
            @GeometricBrownianDriftWithJacobian, ...
            @GeometricBrownianDiffusion, ...
            Tcoarse, ...
            x0, ...
            dWcoarse(i,:), ...
            p);

        XexpFinal(i) = Xexp(end);
        XimpFinal(i) = Ximp(end);

        XexactFinal(i) = x0*exp( ...
            (lambda - 0.5*sigma^2)*Tfinal ...
            + sigma*Wcoarse(i,end));
    end

    errExplicit(j) = mean(abs(XexactFinal-XexpFinal));
    errImplicit(j) = mean(abs(XexactFinal-XimpFinal));
end

fitExplicit = polyfit(log(dtList),log(errExplicit),1);
fitImplicit = polyfit(log(dtList),log(errImplicit),1);

orderExplicit = fitExplicit(1);
orderImplicit = fitImplicit(1);

fprintf('\nEstimated strong convergence orders\n')
fprintf('Explicit-explicit: %.4f\n',orderExplicit)
fprintf('Implicit-explicit: %.4f\n',orderImplicit)

%% Convergence-order plot
figure

loglog(dtList,errExplicit,'o-','LineWidth',1.2)
hold on
loglog(dtList,errImplicit,'s-','LineWidth',1.2)

xlabel('\Deltat')
ylabel('E(|x(T)-x_N|)')
legend( ...
    sprintf('Explicit-explicit, p = %.3f',orderExplicit), ...
    sprintf('Implicit-explicit, p = %.3f',orderImplicit), ...
    'Location','best')
title('Strong convergence order')
grid on

%% Scalar standard Wiener process
function [W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed)

    if nargin == 4
        rng(seed);
    end

    dt = T/N;
    dW = sqrt(dt)*randn(Ns,N);
    W = [zeros(Ns,1) cumsum(dW,2)];
    Tw = 0:dt:T;
end

%% Explicit-explicit Euler-Maruyama SDE solver
function X = SDEeulerExplicitExplicit(ffun,gfun,T,x0,W,varargin)

    N = size(T,2)-1;
    nx = size(x0,1);

    X = zeros(nx,N+1);
    X(:,1) = x0;

    for k = 1:N
        dt = T(k+1)-T(k);
        dW = W(:,k+1)-W(:,k);

        f = feval(ffun,T(k),X(:,k),varargin{:});
        g = feval(gfun,T(k),X(:,k),varargin{:});

        X(:,k+1) = X(:,k) + f*dt + g*dW;
    end
end

%% Newton solver copied from Lecture 03A/03B implicit-explicit structure
function [x,f,J] = SDENewtonSolver(ffun,t,dt,psi,xinit,tol,maxit,varargin)

    I = eye(length(xinit));
    x = xinit;

    [f,J] = feval(ffun,t,x,varargin{:});
    R = x - f*dt - psi;

    it = 1;

    while ((norm(R,'inf') > tol) && (it <= maxit))
        dRdx = I - J*dt;
        mdx = dRdx\R;

        x = x - mdx;

        [f,J] = feval(ffun,t,x,varargin{:});
        R = x - f*dt - psi;

        it = it + 1;
    end

    if norm(R,'inf') > tol
        error('Newton iteration did not converge.')
    end
end

%% Implicit-explicit Euler SDE solver
function X = SDEeulerImplicitExplicit(ffun,gfun,T,x0,dW,varargin)

    tol = 1.0e-8;
    maxit = 10;

    N = size(T,2)-1;
    nx = size(x0,1);

    X = zeros(nx,N+1);
    X(:,1) = x0;

    for k = 1:N
        dt = T(k+1)-T(k);

        f = feval(ffun,T(k),X(:,k),varargin{:});
        g = feval(gfun,T(k),X(:,k),varargin{:});

        psi = X(:,k) + g*dW(:,k);
        xinit = psi + f*dt;

        [X(:,k+1),f,~] = SDENewtonSolver( ...
            ffun, ...
            T(k+1), ...
            dt, ...
            psi, ...
            xinit, ...
            tol, ...
            maxit, ...
            varargin{:});
    end
end

%% Geometric Brownian Motion drift
function f = GeometricBrownianDrift(t,x,p)

    lambda = p(1);
    f = lambda*x;
end

%% Geometric Brownian Motion drift with Jacobian
function [f,J] = GeometricBrownianDriftWithJacobian(t,x,p)

    lambda = p(1);

    f = lambda*x;
    J = lambda;
end

%% Geometric Brownian Motion diffusion
function g = GeometricBrownianDiffusion(t,x,p)

    sigma = p(2);
    g = sigma*x;
end
