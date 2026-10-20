%% Simulation 1: F1,F2 constant

% x = [h1; h2; h3; h4]
% u = [F1; F2]
% p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2]


p = [1.2272; 1.2272; 1.2272; 1.2272; ...
     380.1327; 380.1327; 380.1327; 380.1327; ...
     981; 0.65; 0.55];

% gamma values are minimum-phase values
% a larger fraction of each pump flow goes directly to the lower tanks
% which is what we directly control

t0 = 0;
tf = 20*60;

h10 = 0;
h20 = 0;
h30 = 0;
h40 = 0;

F1 = 300;
F2 = 300;

x0 = [h10; h20; h30; h40];
u = [F1; F2];


% use non-stiff ODE solver

ODEoptions = odeset();
[T,X] = ode45(@QuadrupleTankProcess, ...
               [t0 tf], x0, ODEoptions, u, p);

figure

subplot(2,2,1)
plot(T/60,X(:,1))
xlabel('Time [min]')
ylabel('h_1 [cm]')
grid on

subplot(2,2,2)
plot(T/60,X(:,2))
xlabel('Time [min]')
ylabel('h_2 [cm]')
grid on

subplot(2,2,3)
plot(T/60,X(:,3))
xlabel('Time [min]')
ylabel('h_3 [cm]')
grid on

subplot(2,2,4)
plot(T/60,X(:,4))
xlabel('Time [min]')
ylabel('h_4 [cm]')
grid on

%% PID Setup

F1s = 300;
F2s = 300;

xss = X(end,:)';      % steady state for F1 = F2 = 300

Ts = 5;               % controller sample time [s]
tf_cl = 60*60;        % simulation time [s]

t = 0:Ts:tf_cl;
N = length(t);

% Allocate arrays
x = zeros(4,N);       % states: h1,h2,h3,h4
u = zeros(2,N);       % inputs: F1,F2

% Start closed-loop simulation at nominal equilibrium
x(:,1) = xss;

% Control variables are the first two tank levels
h1ss = xss(1);
h2ss = xss(2);

% Reference starts out at steady state tank levels
r1 = h1ss*ones(1,N);
r2 = h2ss*ones(1,N);

% Change h1 setpoint after 10 minutes
r1(t >= 10*60) = h1ss + 3;

% Change h2 setpoint after 30 minutes
r2(t >= 30*60) = h2ss + 3;

% Bounds on water force
Fmin = 0;
Fmax = 500;

% Initial PID values from slides
I1 = 0;
I2 = 0;

% Previously measured tank levels for use in derivative control
h1_old = x(1,1);
h2_old = x(2,1);

% Values for the two SISO PID controllers
% We control tank 1 with one controller, and tank 2 with another

KP1 = 14;
KI1 = 0.33;
KD1 = 0;

KP2 = 18;
KI2 = 0.39;
KD2 = 0;

%% Closed loop simulation

for k = 1:N-1

    % Current measured levels
    h1 = x(1,k);
    h2 = x(2,k);

    % PID controller for tank 1
    [F1,I1] = SISOPID( ...
        F1s, r1(k), h1, h1_old, I1, ...
        KP1, KI1, KD1, Ts, ...
        Fmin, Fmax);

    % PID controller for tank 2
    [F2,I2] = SISOPID( ...
        F2s, r2(k), h2, h2_old, I2, ...
        KP2, KI2, KD2, Ts, ...
        Fmin, Fmax);

    % Store manipulated inputs
    u(:,k) = [F1; F2];

    % Simulate plant over one controller sample
    [~,Xtemp] = ode15s(@QuadrupleTankProcess, ...
                       [t(k) t(k+1)], ...
                       x(:,k), ...
                       ODEoptions, ...
                       u(:,k), ...
                       p);

    % State at next sampling instant
    x(:,k+1) = Xtemp(end,:)';

    % Previous measurements for derivative term
    h1_old = h1;
    h2_old = h2;

end

%At final plotted time, assume the input is equal to last computed input.
u(:,N) = u(:,N-1);

%% Plot

figure

subplot(2,1,1)
plot(t/60,x(1,:),'LineWidth',1.5)
hold on
plot(t/60,r1,'--','LineWidth',1.2)
xlabel('Time [min]')
ylabel('h_1 [cm]')
legend('h_1','Setpoint')
grid on

subplot(2,1,2)
plot(t/60,x(2,:),'LineWidth',1.5)
hold on
plot(t/60,r2,'--','LineWidth',1.2)
xlabel('Time [min]')
ylabel('h_2 [cm]')
legend('h_2','Setpoint')
grid on

figure

subplot(2,1,1)
plot(t/60,u(1,:),'LineWidth',1.5)
xlabel('Time [min]')
ylabel('F_1 [cm^3/s]')
grid on

subplot(2,1,2)
plot(t/60,u(2,:),'LineWidth',1.5)
xlabel('Time [min]')
ylabel('F_2 [cm^3/s]')
grid on

%% Stress testing


F1s = 300;
F2s = 300;

xss = X(end,:)';      % steady state for F1 = F2 = 300

Ts = 5;               % controller sample time [s]
tf_cl = 60*60;        % simulation time [s]

t = 0:Ts:tf_cl;
N = length(t);

% Allocate arrays
x = zeros(4,N);       % states: h1,h2,h3,h4
u = zeros(2,N);       % inputs: F1,F2

% Start closed-loop simulation at nominal equilibrium
x(:,1) = xss;

% Control variables are the first two tank levels
h1ss = xss(1);
h2ss = xss(2);

% Reference starts at steady state
r1 = h1ss*ones(1,N);
r2 = h2ss*ones(1,N);

% Frequent setpoint changes for stress testing

% Tank 1
r1(t >=  5*60 & t < 10*60) = h1ss + 3;
r1(t >= 10*60 & t < 15*60) = h1ss - 2;
r1(t >= 15*60 & t < 20*60) = h1ss + 5;
r1(t >= 20*60 & t < 25*60) = h1ss;
r1(t >= 25*60 & t < 30*60) = h1ss - 4;
r1(t >= 30*60 & t < 35*60) = h1ss + 2;
r1(t >= 35*60 & t < 40*60) = h1ss + 6;
r1(t >= 40*60 & t < 45*60) = h1ss - 3;
r1(t >= 45*60 & t < 50*60) = h1ss + 1;
r1(t >= 50*60)              = h1ss;

% Tank 2
r2(t >=  5*60 & t < 10*60) = h2ss - 2;
r2(t >= 10*60 & t < 15*60) = h2ss + 4;
r2(t >= 15*60 & t < 20*60) = h2ss - 3;
r2(t >= 20*60 & t < 25*60) = h2ss + 5;
r2(t >= 25*60 & t < 30*60) = h2ss;
r2(t >= 30*60 & t < 35*60) = h2ss - 4;
r2(t >= 35*60 & t < 40*60) = h2ss + 3;
r2(t >= 40*60 & t < 45*60) = h2ss + 6;
r2(t >= 45*60 & t < 50*60) = h2ss - 2;
r2(t >= 50*60)              = h2ss;

% Bounds on water force
Fmin = 0;
Fmax = 500;

% Initial PID values from slides
I1 = 0;
I2 = 0;

% Previously measured tank levels for use in derivative control
h1_old = x(1,1);
h2_old = x(2,1);

% Values for the two SISO PID controllers
% We control tank 1 with one controller, and tank 2 with another

KP1 = 14;
KI1 = 0.33;
KD1 = 0;

KP2 = 18;
KI2 = 0.39;
KD2 = 0;

for k = 1:N-1

    % Current measured levels
    h1 = x(1,k);
    h2 = x(2,k);

    % PID controller for tank 1
    [F1,I1] = SISOPID( ...
        F1s, r1(k), h1, h1_old, I1, ...
        KP1, KI1, KD1, Ts, ...
        Fmin, Fmax);

    % PID controller for tank 2
    [F2,I2] = SISOPID( ...
        F2s, r2(k), h2, h2_old, I2, ...
        KP2, KI2, KD2, Ts, ...
        Fmin, Fmax);

    % Store manipulated inputs
    u(:,k) = [F1; F2];

    % Simulate plant over one controller sample
    [~,Xtemp] = ode15s(@QuadrupleTankProcess, ...
                       [t(k) t(k+1)], ...
                       x(:,k), ...
                       ODEoptions, ...
                       u(:,k), ...
                       p);

    % State at next sampling instant
    x(:,k+1) = Xtemp(end,:)';

    % Previous measurements for derivative term
    h1_old = h1;
    h2_old = h2;

end

%At final plotted time, assume the input is equal to last computed input.
u(:,N) = u(:,N-1);

figure

subplot(2,1,1)
plot(t/60,x(1,:),'LineWidth',1.5)
hold on
plot(t/60,r1,'--','LineWidth',1.2)
xlabel('Time [min]')
ylabel('h_1 [cm]')
legend('h_1','Setpoint')
grid on

subplot(2,1,2)
plot(t/60,x(2,:),'LineWidth',1.5)
hold on
plot(t/60,r2,'--','LineWidth',1.2)
xlabel('Time [min]')
ylabel('h_2 [cm]')
legend('h_2','Setpoint')
grid on

figure

subplot(2,1,1)
plot(t/60,u(1,:),'LineWidth',1.5)
xlabel('Time [min]')
ylabel('F_1 [cm^3/s]')
grid on

subplot(2,1,2)
plot(t/60,u(2,:),'LineWidth',1.5)
xlabel('Time [min]')
ylabel('F_2 [cm^3/s]')
grid on