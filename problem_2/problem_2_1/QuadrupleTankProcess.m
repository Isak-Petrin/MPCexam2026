function xdot = QuadrupleTankProcess(~,x,u,d,p)
% Modified Four Tank System: xdot = f(x,u,d,p)
% x = [m1; m2; m3; m4] [g]       masses in the four tanks
% u = [F1; F2]         [cm^3/s]  manipulated pump flows
% d = [F3; F4]         [cm^3/s]  unmeasured disturbance flows into tank 3 and 4

m = x;
F = u;

a     = p(1:4,1);     % [cm^2]   pipe cross sectional areas
A     = p(5:8,1);     % [cm^2]   tank cross sectional areas
gamma = p(9:10,1);    % [-]      valve constants
g     = p(11,1);      % [cm/s^2] acceleration of gravity
rho   = p(12,1);      % [g/cm^3] density of water

% Inflows
qin = zeros(4,1);
qin(1) = gamma(1)*F(1);
qin(2) = gamma(2)*F(2);
qin(3) = (1-gamma(2))*F(2) + d(1);   % F3 enters tank 3
qin(4) = (1-gamma(1))*F(1) + d(2);   % F4 enters tank 4

% Levels (guard against negative mass during noisy simulations)
h = max(m./(rho*A), 0);

% Outflows
qout = a.*sqrt(2*g*h);

% Differential equations (mass balances)
xdot = zeros(4,1);
xdot(1) = rho*(qin(1) + qout(3) - qout(1));
xdot(2) = rho*(qin(2) + qout(4) - qout(2));
xdot(3) = rho*(qin(3) - qout(3));
xdot(4) = rho*(qin(4) - qout(4));
end
