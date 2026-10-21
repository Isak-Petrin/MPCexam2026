%% Exercise 3
% Fixed-step explicit Euler
% Fixed-step implicit Euler
% Adaptive-step explicit Euler
% Adaptive-step implicit Euler
% Test all four solvers on the Van der Pol oscillator

clear
clc
close all

%% Van der Pol setup from Lecture 03B
mu = 3;
x0 = [0.5; 0.5];
p = [mu];

tf = 5*mu;

%% Fixed step size
dtFixed = 0.01;
Tfixed = 0:dtFixed:tf;

XexplicitFixed = ODEsolverExplicitEulerFixed( ...
    @VanderpolDrift,Tfixed,x0,p);

XimplicitFixed = ODEsolverImplicitEulerFixed( ...
    @VanderpolDrift,Tfixed,x0,p);

%% Adaptive step size
dt0 = 0.01;
dtMin = 1.0e-6;
dtMax = 0.25;

atol = 1.0e-6;
rtol = 1.0e-4;

[TexplicitAdaptive,XexplicitAdaptive] = ODEsolverExplicitEulerAdaptive( ...
    @VanderpolDrift,[0 tf],x0,dt0,dtMin,dtMax,atol,rtol,p);

[TimplicitAdaptive,XimplicitAdaptive] = ODEsolverImplicitEulerAdaptive( ...
    @VanderpolDrift,[0 tf],x0,dt0,dtMin,dtMax,atol,rtol,p);

%% Reference solution for comparison
odeOptions = odeset('RelTol',1.0e-10,'AbsTol',1.0e-12);

[Tref,Xref] = ode45( ...
    @(t,x) VanderpolDrift(t,x,p), ...
    [0 tf], ...
    x0, ...
    odeOptions);

Xref = Xref';

%% Time-series plots
figure

subplot(2,1,1)
plot(Tfixed,XexplicitFixed(1,:),'LineWidth',1.0)
hold on
plot(Tfixed,XimplicitFixed(1,:),'LineWidth',1.0)
plot(TexplicitAdaptive,XexplicitAdaptive(1,:),'--','LineWidth',1.0)
plot(TimplicitAdaptive,XimplicitAdaptive(1,:),'--','LineWidth',1.0)
plot(Tref,Xref(1,:),'k','LineWidth',1.5)
xlabel('t')
ylabel('x_1')
legend('Explicit fixed','Implicit fixed', ...
       'Explicit adaptive','Implicit adaptive', ...
       'ode45 reference','Location','best')
grid on

subplot(2,1,2)
plot(Tfixed,XexplicitFixed(2,:),'LineWidth',1.0)
hold on
plot(Tfixed,XimplicitFixed(2,:),'LineWidth',1.0)
plot(TexplicitAdaptive,XexplicitAdaptive(2,:),'--','LineWidth',1.0)
plot(TimplicitAdaptive,XimplicitAdaptive(2,:),'--','LineWidth',1.0)
plot(Tref,Xref(2,:),'k','LineWidth',1.5)
xlabel('t')
ylabel('x_2')
legend('Explicit fixed','Implicit fixed', ...
       'Explicit adaptive','Implicit adaptive', ...
       'ode45 reference','Location','best')
grid on

%% Phase-plane plot
figure

plot(XexplicitFixed(1,:),XexplicitFixed(2,:),'LineWidth',1.0)
hold on
plot(XimplicitFixed(1,:),XimplicitFixed(2,:),'LineWidth',1.0)
plot(XexplicitAdaptive(1,:),XexplicitAdaptive(2,:),'--','LineWidth',1.0)
plot(XimplicitAdaptive(1,:),XimplicitAdaptive(2,:),'--','LineWidth',1.0)
plot(Xref(1,:),Xref(2,:),'k','LineWidth',1.5)

xlabel('x_1')
ylabel('x_2')
legend('Explicit fixed','Implicit fixed', ...
       'Explicit adaptive','Implicit adaptive', ...
       'ode45 reference','Location','best')
grid on

%% Adaptive step sizes
figure

subplot(2,1,1)
stairs(TexplicitAdaptive(1:end-1),diff(TexplicitAdaptive),'LineWidth',1.0)
xlabel('t')
ylabel('\Deltat')
title('Adaptive explicit Euler step size')
grid on

subplot(2,1,2)
stairs(TimplicitAdaptive(1:end-1),diff(TimplicitAdaptive),'LineWidth',1.0)
xlabel('t')
ylabel('\Deltat')
title('Adaptive implicit Euler step size')
grid on

%% Number of accepted steps
fprintf('Fixed-step explicit Euler steps:   %d\n',length(Tfixed)-1)
fprintf('Fixed-step implicit Euler steps:   %d\n',length(Tfixed)-1)
fprintf('Adaptive explicit Euler steps:     %d\n',length(TexplicitAdaptive)-1)
fprintf('Adaptive implicit Euler steps:     %d\n',length(TimplicitAdaptive)-1)

%% Fixed-step explicit Euler
function X = ODEsolverExplicitEulerFixed(ffun,T,x0,varargin)

    N = length(T);
    nx = length(x0);

    X = zeros(nx,N);
    X(:,1) = x0;

    for k = 1:N-1
        f = feval(ffun,T(k),X(:,k),varargin{:});
        dt = T(k+1)-T(k);

        X(:,k+1) = X(:,k) + f*dt;
    end
end

%% Fixed-step implicit Euler
function X = ODEsolverImplicitEulerFixed(ffun,T,x0,varargin)

    tol = 1.0e-10;
    maxit = 100;

    N = length(T);
    nx = length(x0);

    X = zeros(nx,N);
    X(:,1) = x0;

    for k = 1:N-1
        dt = T(k+1)-T(k);

        [f,~] = feval(ffun,T(k),X(:,k),varargin{:});
        xinit = X(:,k) + f*dt;

        X(:,k+1) = ODENewtonSolver( ...
            ffun,T(k+1),dt,X(:,k),xinit,tol,maxit,varargin{:});
    end
end

%% Newton solver copied from the implicit solver structure in Lecture 03B
function [x,f,J] = ODENewtonSolver(ffun,t,dt,psi,xinit,tol,maxit,varargin)

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

%% Adaptive-step explicit Euler
function [T,X] = ODEsolverExplicitEulerAdaptive( ...
    ffun,tspan,x0,dt0,dtMin,dtMax,atol,rtol,varargin)

    t = tspan(1);
    tf = tspan(2);
    x = x0;
    dt = dt0;

    T = t;
    X = x;

    safety = 0.9;
    facMin = 0.2;
    facMax = 5.0;

    while t < tf

        dt = min(dt,tf-t);

        f = feval(ffun,t,x,varargin{:});
        xFull = x + dt*f;

        dtHalf = dt/2;

        xHalf = x + dtHalf*f;

        fHalf = feval(ffun,t+dtHalf,xHalf,varargin{:});
        xTwoHalf = xHalf + dtHalf*fHalf;

        err = norm(xTwoHalf-xFull,'inf');

        %apparently this one is standard procedure
        scale = atol + rtol*max(norm(x,'inf'),norm(xTwoHalf,'inf'));

        if err <= scale || dt <= dtMin
            t = t + dt;
            x = xTwoHalf;

            T(end+1) = t;
            X(:,end+1) = x;
        end

        if err == 0
            fac = facMax;
        else
            fac = safety*(scale/err)^(1/2);
            fac = min(facMax,max(facMin,fac));
        end

        dt = min(dtMax,max(dtMin,fac*dt));
    end
end

%% Adaptive-step implicit Euler
function [T,X] = ODEsolverImplicitEulerAdaptive( ...
    ffun,tspan,x0,dt0,dtMin,dtMax,atol,rtol,varargin)

    t = tspan(1);
    tf = tspan(2);
    x = x0;
    dt = dt0;

    T = t;
    X = x;

    newtonTol = 1.0e-10;
    maxit = 100;

    safety = 0.9;
    facMin = 0.2;
    facMax = 5.0;

    while t < tf

        dt = min(dt,tf-t);

        xFull = ImplicitEulerStep( ...
            ffun,t,x,dt,newtonTol,maxit,varargin{:});

        dtHalf = dt/2;

        xHalf = ImplicitEulerStep( ...
            ffun,t,x,dtHalf,newtonTol,maxit,varargin{:});

        xTwoHalf = ImplicitEulerStep( ...
            ffun,t+dtHalf,xHalf,dtHalf,newtonTol,maxit,varargin{:});

        err = norm(xTwoHalf-xFull,'inf');

        %apparently this one is standard procedure
        scale = atol + rtol*max(norm(x,'inf'),norm(xTwoHalf,'inf'));

        if err <= scale || dt <= dtMin
            t = t + dt;
            x = xTwoHalf;

            T(end+1) = t;
            X(:,end+1) = x;
        end

        if err == 0
            fac = facMax;
        else
            fac = safety*(scale/err)^(1/2);
            fac = min(facMax,max(facMin,fac));
        end

        dt = min(dtMax,max(dtMin,fac*dt));
    end
end

%% One implicit Euler step
function xNext = ImplicitEulerStep(ffun,t,x,dt,tol,maxit,varargin)

    [f,~] = feval(ffun,t,x,varargin{:});
    xinit = x + f*dt;

    xNext = ODENewtonSolver( ...
        ffun,t+dt,dt,x,xinit,tol,maxit,varargin{:});
end

%% Van der Pol drift copied from Lecture 03B
function [f,J] = VanderpolDrift(t,x,p)

    mu = p(1);
    tmp = mu*(1.0-x(1)*x(1));

    f = zeros(2,1);
    f(1,1) = x(2);
    f(2,1) = tmp*x(2)-x(1);

    if nargout > 1
        J = [0 1; -2*mu*x(1)*x(2)-1.0 tmp];
    end
end
