% Functions in slide order:
% ScalarStdWienerProcess -> ScalarSampleMeanStdVar -> SDEeulerExplicitExplicit -> LangevinDrift -> LangevinDiffusion -> GeometricBrownianDrift -> GeometricBrownianDiffusion -> SDENewtonSolver -> SDEeulerImplicitExplicit -> StdWienerProcess -> SDEsolverExplicitExplicit -> SDENewtonSolver (repeated in Lecture 03B) -> SDEsolverImplicitExplicit -> VanderpolDrift -> VanderPolDiffusion1 -> VanderPolDiffusion2
%
% This file follows the order in which MATLAB implementations appear in Lecture 03A, then Lecture 03B.
% Repeated implementations are shown again where they reappear in the slides.
% Because MATLAB does not allow duplicate local function names in one executable file, this file is intended as a slide-ordered reference rather than a single run-all script.

%% Lecture 03A - slide 9: Realizations of Standard Brownian Motion
Ns = 10;
T = 10;
N = 1000;
seed = 100;

rng(seed);
dt = T/N;
dW = sqrt(dt)*randn(Ns,N);
W = cumsum(dW,2);

%% Lecture 03A - slide 10: Scalar Standard Wiener Process
function [W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed)

    if nargin == 4
        rng(seed);
    end

    dt = T/N;
    dW = sqrt(dt)*randn(Ns,N);
    W = [zeros(Ns,1) cumsum(dW,2)];
    Tw = 0:dt:T;
end

[W,Tw,dW] = ScalarStdWienerProcess(10,1000,10,100);
plot(Tw,W,'linewidth',2);

%% Lecture 03A - slide 14: Scalar Sample Mean and Standard Deviation
function [xmean,s,xmeanp2s,xmeanm2s] = ScalarSampleMeanStdVar(x)

    xmean = mean(x,1);
    s = std(x,0,1);
    xmeanp2s = xmean + 2*s;
    xmeanm2s = xmean - 2*s;
end

[W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed);
[Wmean,sW,Wmeanp2sW,Wmeanm2sW] = ScalarSampleMeanStdVar(W);

%% Lecture 03A - slide 28: Langevin direct MATLAB implementation
T = 10;
N = 1000;
Ns = 10;
seed = 100;

x0 = 10;
lambda = -0.5;
sigma = 1.0;

[W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed);

X = zeros(size(W));

for i = 1:Ns
    X(i,1) = x0;

    for k = 1:N
        dt = Tw(k+1)-Tw(k);

        X(i,k+1) = X(i,k) ...
            + lambda*X(i,k)*dt ...
            + sigma*dW(i,k);
    end
end

%% Lecture 03A - slide 32: Explicit-explicit SDE solver
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

%% Lecture 03A - slide 33: Langevin driver, drift and diffusion
T = 10;
N = 1000;
Ns = 10;
seed = 100;

x0 = 10;
lambda = -0.5;
sigma = 1.0;
p = [lambda; sigma];

[W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed);

X = zeros(size(W));

for i = 1:Ns
    X(i,:) = SDEeulerExplicitExplicit( ...
        @LangevinDrift, ...
        @LangevinDiffusion, ...
        Tw, ...
        x0, ...
        W(i,:), ...
        p);
end

function f = LangevinDrift(t,x,p)

    lambda = p(1);
    f = lambda*x;
end

function g = LangevinDiffusion(t,x,p)

    sigma = p(2);
    g = sigma;
end

%% Lecture 03A - slide 35: Geometric Brownian Motion implementation
T = 10;
N = 1000;
Ns = 10;
seed = 100;

x0 = 1;
lambda = 0.15;
sigma = 0.15;
p = [lambda; sigma];

[W,Tw,dW] = ScalarStdWienerProcess(T,N,Ns,seed);

X = zeros(size(W));

for i = 1:Ns
    X(i,:) = SDEeulerExplicitExplicit( ...
        @GeometricBrownianDrift, ...
        @GeometricBrownianDiffusion, ...
        Tw, ...
        x0, ...
        W(i,:), ...
        p);
end

function f = GeometricBrownianDrift(t,x,p)

    lambda = p(1);
    f = lambda*x;
end

function g = GeometricBrownianDiffusion(t,x,p)

    sigma = p(2);
    g = sigma*x;
end

%% Lecture 03A - slide 40: Newton method in implicit-explicit SDE solver
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
end

%% Lecture 03A - slide 41: Implicit-explicit SDE solver
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
            T(:,k+1), ...
            dt, ...
            psi, ...
            xinit, ...
            tol, ...
            maxit, ...
            varargin{:});
    end
end

%% Lecture 03B - slide 6: Multivariate Standard Wiener Process
function [W,Tw,dW] = StdWienerProcess(T,N,nW,Ns,seed)

    if nargin == 4
        rng(seed);
    end

    dt = T/N;

    dW = sqrt(dt)*randn(nW,N,Ns);
    W = [zeros(nW,1,Ns) cumsum(dW,2)];

    Tw = 0:dt:T;
end

%% Lecture 03B - slide 7: Explicit-explicit SDE solver
function X = SDEsolverExplicitExplicit(ffun,gfun,T,x0,W,varargin)

    N = length(T);
    nx = length(x0);

    X = zeros(nx,N);
    X(:,1) = x0;

    for k = 1:N-1
        f = feval(ffun,T(k),X(:,k),varargin{:});
        g = feval(gfun,T(k),X(:,k),varargin{:});

        dt = T(k+1)-T(k);
        dW = W(:,k+1)-W(:,k);

        psi = X(:,k) + g*dW;

        X(:,k+1) = psi + f*dt;
    end
end

%% Lecture 03B - slide 11: Newton method in implicit-explicit SDE solver
function [x,f,J] = SDENewtonSolver_L03B(ffun,t,dt,psi,xinit,tol,maxit,varargin)

    I = eye(length(xinit));
    x = xinit;

    [f,J] = feval(ffun,t,x,varargin{:});
    R = x - f*dt - psi;

    it = 1;-

    while ((norm(R,'inf') > tol) && (it <= maxit))
        dRdx = I - J*dt;
        mdx = dRdx\R;

        x = x - mdx;

        [f,J] = feval(ffun,t,x,varargin{:});
        R = x - f*dt - psi;

        it = it + 1;
    end
end

%% Lecture 03B - slide 12: Implicit-explicit SDE solver
function X = SDEsolverImplicitExplicit(ffun,gfun,T,x0,W,varargin)

    tol = 1.0e-8;
    maxit = 100;

    N = length(T);
    nx = length(x0);

    X = zeros(nx,N);
    X(:,1) = x0;

    k = 1;
    [f,~] = feval(ffun,T(k),X(:,k),varargin{:});

    for k = 1:N-1
        g = feval(gfun,T(k),X(:,k),varargin{:});

        dt = T(k+1)-T(k);
        dW = W(:,k+1)-W(:,k);

        psi = X(:,k) + g*dW;
        xinit = psi + f*dt;

        [X(:,k+1),f,~] = SDENewtonSolver_L03B( ...
            ffun, ...
            T(:,k+1), ...
            dt, ...
            psi, ...
            xinit, ...
            tol, ...
            maxit, ...
            varargin{:});
    end
end

%% Lecture 03B - slide 14: Van der Pol drift and diffusion functions
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

function g = VanderPolDiffusion1(t,x,p)

    sigma = p(2);
    g = [0.0; sigma];
end

function g = VanderPolDiffusion2(t,x,p)

    sigma = p(2);
    g = [0.0; sigma*(1.0+x(1)*x(1))];
end

%% Lecture 03B - slide 17: Van der Pol state-independent diffusion driver
mu = 3;
sigma = 1.0;

x0 = [0.5; 0.5];
p = [mu; sigma];

tf = 5*mu;
nw = 1;
N = 1000;
Ns = 5;
seed = 100;

[W,T,~] = StdWienerProcess(tf,N,nw,Ns,seed);

X = zeros(length(x0),N+1,Ns);

for i = 1:Ns
    X(:,:,i) = SDEsolverExplicitExplicit( ...
        @VanderpolDrift, ...
        @VanderPolDiffusion1, ...
        T, ...
        x0, ...
        W(:,:,i), ...
        p);
end

Xd = SDEsolverExplicitExplicit( ...
    @VanderpolDrift, ...
    @VanderPolDiffusion1, ...
    T, ...
    x0, ...
    W(:,:,i), ...
    [mu; 0.0]);

%% Lecture 03B - slide 20: Van der Pol state-dependent diffusion driver
mu = 3;
sigma = 0.5;

x0 = [0.5; 0.5];
p = [mu; sigma];

tf = 5*mu;
nw = 1;
N = 5000;
Ns = 5;
seed = 100;

[W,T,~] = StdWienerProcess(tf,N,nw,Ns,seed);

X = zeros(length(x0),N+1,Ns);

for i = 1:Ns
    X(:,:,i) = SDEsolverExplicitExplicit( ...
        @VanderpolDrift, ...
        @VanderPolDiffusion2, ...
        T, ...
        x0, ...
        W(:,:,i), ...
        p);
end

Xd = SDEsolverExplicitExplicit( ...
    @VanderpolDrift, ...
    @VanderPolDiffusion2, ...
    T, ...
    x0, ...
    W(:,:,i), ...
    [mu; 0.0]);
