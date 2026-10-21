%% Deterministic simulation  % Slides 18-25

a1 = 1.2272;   % [cm^2]
a2 = 1.2272;   % [cm^2]
a3 = 1.2272;   % [cm^2]
a4 = 1.2272;   % [cm^2]

A1 = 380.1327; % [cm^2]
A2 = 380.1327; % [cm^2]
A3 = 380.1327; % [cm^2]
A4 = 380.1327; % [cm^2]

gamma1 = 0.45;
gamma2 = 0.40;

g   = 981;     % [cm/s^2]
rho = 1.00;    % [g/cm^3]

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];


% time grid: t_k = 60k seconds, k = 0,...20
Ts_min = 1;            % [min]
Ts     = 60;           % [s]
tmin   = 0:Ts_min:20;  % [min]
t      = tmin*60;      % [s]
N      = length(t);


nx = 4; % 4 states
ny = 4; % 4 tank levels measured
nz = 2; % 2 tank levels controlled

N = length(t); % number of sampling instants


x = zeros(nx,N); % each column is one sampled state
y = zeros(ny,N); % stores sensor measurements at each instant
z = zeros(nz,N); % stores controlled outputs 


% intermediate values created by the ODE solver
% if we sample at t=10 and t=20 the ODE solver perhaps approximates values
% at t = 15 too.

X = zeros(0,nx);
T = zeros(0,1);

x0 = [0; 0; 0; 0]; %initial mass
x(:,1) = x0; % initial condition on states



F1 = [330 235 295 303 242 292 300 243 325 324 ...
      330 357 278 274 202 301 320 334 276 308];

F2 = [285 301 295 313 314 316 289 278 302 293 ...
      304 296 307 277 319 305 300 306 307 297];

u = [F1; F2];

for k = 1:N-1

    % Sensor values (all 4 levels)
    y(:,k) = QuadrupleTankSensor(x(:,k),p);

    % Controlled outputs (lower tank levels)
    z(:,k) = QuadrupleTankOutput(x(:,k),p);

    % Simulate one sampling interval
    [Tk,Xk] = ode15s(@QuadrupleTankProcess, ...
                     [t(k) t(k+1)], ...
                     x(:,k), ...
                     [], ...
                     u(:,k), ...
                     p);

    % State at next sample
    x(:,k+1) = Xk(end,:)';

    % Store detailed trajectory
    T = [T; Tk];
    X = [X; Xk];
end


% Final sample
k = N;

y(:,k) = QuadrupleTankSensor(x(:,k),p);
z(:,k) = QuadrupleTankOutput(x(:,k),p);


figure;
subplot(1,2,1)
stairs(tmin,[F1 F1(end)],'LineWidth',1.5)
xlabel('time [min]')
ylabel('F_1 [cm^3/s]')
title('Pump Flow F_1')
grid on
xlim([0 20])

subplot(1,2,2)
stairs(tmin,[F2 F2(end)],'LineWidth',1.5)
xlabel('time [min]')
ylabel('F_2 [cm^3/s]')
title('Pump Flow F_2')
grid on
xlim([0 20])

figure;
subplot(2,2,1)
plot(tmin,x(3,:)/1000,'o-','LineWidth',1.2)
xlabel('time [min]')
ylabel('m_3 [kg]')
title('Tank 3')
grid on
xlim([0 20])

subplot(2,2,2)
plot(tmin,x(4,:)/1000,'o-','LineWidth',1.2)
xlabel('time [min]')
ylabel('m_4 [kg]')
title('Tank 4')
grid on
xlim([0 20])

subplot(2,2,3)
plot(tmin,x(1,:)/1000,'o-','LineWidth',1.2)
xlabel('time [min]')
ylabel('m_1 [kg]')
title('Tank 1')
grid on
xlim([0 20])

subplot(2,2,4)
plot(tmin,x(2,:)/1000,'o-','LineWidth',1.2)
xlabel('time [min]')
ylabel('m_2 [kg]')
title('Tank 2')
grid on
xlim([0 20])

figure;
subplot(2,2,1)
plot(tmin,y(3,:),'o-','LineWidth',1.2)
xlabel('time [min]')
ylabel('h_3 [cm]')
title('Tank 3')
grid on
xlim([0 20])

subplot(2,2,2)
plot(tmin,y(4,:),'o-','LineWidth',1.2)
xlabel('time [min]')
ylabel('h_4 [cm]')
title('Tank 4')
grid on
xlim([0 20])

subplot(2,2,3)
plot(tmin,y(1,:),'o-','LineWidth',1.2)
xlabel('time [min]')
ylabel('h_1 [cm]')
title('Tank 1')
grid on
xlim([0 20])

subplot(2,2,4)
plot(tmin,y(2,:),'o-','LineWidth',1.2)
xlabel('time [min]')
ylabel('h_2 [cm]')
title('Tank 2')
grid on
xlim([0 20])

%% Naive stochastic simulation  % Slides 26-34

% sample every 10 seconds
t0 = 0;
tf = 20*60;
Ts = 10;
t = t0:Ts:tf;
tmin = t/60;
N = length(t);


% constant forces
F1 = 300;      
F2 = 300;       
u = [repmat(F1,1,N);
     repmat(F2,1,N)];


rng(100);

% process noise
Q = [20^2  0;
     0     40^2];
Lq = chol(Q,'lower');
w = Lq*randn(2,N);

% measuremenet noise
R = eye(4);
Lr = chol(R,'lower');
v = Lr*randn(4,N);


for k = 1:N-1

    % noise in sensor
    y(:,k) = QuadrupleTankSensor(x(:,k),p) + v(:,k);

   
    z(:,k) = QuadrupleTankOutput(x(:,k),p);

    
    [Tk,Xk] = ode15s(@QuadrupleTankProcess, ...
                 [t(k) t(k+1)], ...
                 x(:,k), ...
                 [], ...
                 u(:,k) + w(:,k), ... % noise in process
                 p);

    % State at next sample
    x(:,k+1) = Xk(end,:)';

    % Store detailed trajectory
    T = [T; Tk];
    X = [X; Xk];
end

% Final sample
k = N;

y(:,k) = QuadrupleTankSensor(x(:,k),p) + v(:,k);
z(:,k) = QuadrupleTankOutput(x(:,k),p);



Factual = u + w;

figure;

subplot(2,2,1)
plot(tmin,z(1,:),'LineWidth',1.2)
xlabel('time [min]')
ylabel('h_1 [cm]')
title('Tank 1')
grid on

subplot(2,2,2)
plot(tmin,z(2,:),'LineWidth',1.2)
xlabel('time [min]')
ylabel('h_2 [cm]')
title('Tank 2')
grid on

subplot(2,2,3)
stairs(tmin,Factual(1,:),'LineWidth',1.2)
xlabel('time [min]')
ylabel('F_1 [cm^3/s]')
grid on

subplot(2,2,4)
stairs(tmin,Factual(2,:),'LineWidth',1.2)
xlabel('time [min]')
ylabel('F_2 [cm^3/s]')
grid on


ytrue = zeros(4,N);

for k = 1:N
    ytrue(:,k) = QuadrupleTankSensor(x(:,k),p);
end

figure;

subplot(2,2,1)
plot(tmin,y(3,:),'-',tmin,ytrue(3,:),'--','LineWidth',1.2)
xlabel('time [min]')
ylabel('h_3 [cm]')
title('Tank 3')
legend('Noisy measurement','True level','Location','best')
grid on

subplot(2,2,2)
plot(tmin,y(4,:),'-',tmin,ytrue(4,:),'--','LineWidth',1.2)
xlabel('time [min]')
ylabel('h_4 [cm]')
title('Tank 4')
legend('Noisy measurement','True level','Location','best')
grid on

subplot(2,2,3)
plot(tmin,y(1,:),'-',tmin,ytrue(1,:),'--','LineWidth',1.2)
xlabel('time [min]')
ylabel('h_1 [cm]')
title('Tank 1')
legend('Noisy measurement','True level','Location','best')
grid on

subplot(2,2,4)
plot(tmin,y(2,:),'-',tmin,ytrue(2,:),'--','LineWidth',1.2)
xlabel('time [min]')
ylabel('h_2 [cm]')
title('Tank 2')
legend('Noisy measurement','True level','Location','best')
grid on


%% Steady-state  % Slide 35

us = [300;300]; %operating point is F1=F2=300
xs0 = 5000*ones(4,1); %initial guess for solver

xs = fsolve(@(x) QuadrupleTankProcess(0,x,us,p), xs0); %solves f(x,u)=0


ys = QuadrupleTankSensor(xs,p);
zs = QuadrupleTankOutput(xs,p);

disp('us =')
disp(us)

disp('xs =')
disp(xs)

disp('ys =')
disp(ys)

disp('zs =')
disp(zs)


%% Linearization using model - code from slide  % Slides 39-42
% --------------------------------------------------------------
% Parameters
% --------------------------------------------------------------
ap = [a1; a2; a3; a4];       % [cm2] Pipe cross sectional areas
At = [A1; A2; A3; A4];       % [cm2] Tank cross sectional areas
gam = [gamma1; gamma2];       % [-] Valve constants
g = 981;                      % [cm/s2] Acceleration of gravity
rho = 1.00;                   % [g/cm3] Density of water

p = [ap; At; gam; g; rho];


hs = ys;

T = (At./ap).*sqrt(2*hs/g);

A = [-1/T(1)   0         1/T(3)    0;
      0        -1/T(2)   0         1/T(4);
      0         0       -1/T(3)    0;
      0         0        0        -1/T(4)];

B = [rho*gam(1)         0;
     0                  rho*gam(2);
     0                  rho*(1-gam(2));
     rho*(1-gam(1))     0];

C = diag(1./(rho*At));

Cz = C(1:2,:);


%% get 2x2 transfer function matrix  % Slides 43-49
D = zeros(2,2);
clear tf
sys_ss = ss(A,B,Cz,D);
G = tf(sys_ss)
G

%% discrete linear system, numerically  % Slides 50-54 (especially 52-53)

Ts = 10;   % [s] stochastic simulation sampling time

nx = size(A,1);
nu = size(B,2);

% Slide formula:
% [Ad Bd;
%   0  I] = exp([A B;
%                0 0]*Ts)

M = expm([A B;
          zeros(nu,nx) zeros(nu,nu)] * Ts);

Ad = M(1:nx,1:nx);
Bd = M(1:nx,nx+1:nx+nu);

Cd = Cz;
Dd = zeros(size(Cz,1),nu);

disp('Ad =')
disp(Ad)

disp('Bd =')
disp(Bd)

disp('Cd =')
disp(Cd)

disp('Dd =')
disp(Dd)

sys_d = ss(Ad,Bd,Cd,Dd,Ts); %creates state-space object
 
%% direct Matlab command  % Slides 52-53

sys_c = ss(A,B,Cz,zeros(2,2));

Ts = 10;

sys_d = c2d(sys_c,Ts,'zoh');

Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;

sys_d = ss(Ad,Bd,Cd,Dd,Ts);


%% discrete stochastic model  % Slide 70

clc
clear G

sigma = diag(0.10*us);
G = B*sigma;

% x(k+1) = Abar*x(k) + Bbar*u(k) + Gbar*w(k)
% w(k) ~ N(0,Qbar)

% -------------------------------------------------------------
% Abar and Bbar
%
% [Abar Bbar;
%    0    I ] = exp([A B;
%                    0 0] * Ts)
% -------------------------------------------------------------

PhiAB = expm([A B;
              zeros(nu,nx) zeros(nu,nu)] * Ts);

Abar = PhiAB(1:nx,1:nx);
Bbar = PhiAB(1:nx,nx+1:nx+nu);


% -------------------------------------------------------------
% Gbar and Qbar
%
% Gbar = I
%
% Qbar = integral_0^Ts
%        exp(A*tau)*G*G'*exp(A'*tau) dtau
%
% Slide 70 computes Qbar from:
%
% [Phi11 Phi12;
%    0    Phi22]
%
% = exp([-A   G*G';
%         0    A'] * Ts)
%
% Qbar = Phi22' * Phi12
% -------------------------------------------------------------

Gbar = eye(nx);

PhiQ = expm([-A            G*G';
              zeros(nx,nx) A'] * Ts);

Phi12 = PhiQ(1:nx,nx+1:2*nx);
Phi22 = PhiQ(nx+1:2*nx,nx+1:2*nx);

Qbar = Phi22' * Phi12;

% Remove tiny numerical asymmetry caused by floating-point arithmetic
Qbar = (Qbar + Qbar')/2;


%% display Slide 70 matrices

disp('Abar =')
disp(Abar)

disp('Bbar =')
disp(Bbar)

disp('Gbar =')
disp(Gbar)

disp('Qbar =')
disp(Qbar)