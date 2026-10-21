function xdot = QuadrupleTankProcess(~,x,u,p)
% x = [m1; m2; m3; m4] [g]
% u = [F1; F2] [cm^3/s]

m = x;
F = u;

a     = p(1:4,1);     % [cm^2]
A     = p(5:8,1);     % [cm^2]
gamma = p(9:10,1);    % [-]
g     = p(11,1);      % [cm/s^2]
rho   = p(12,1);      % [g/cm^3]

% Inflows
qin = zeros(4,1);
qin(1) = gamma(1)*F(1);
qin(2) = gamma(2)*F(2);
qin(3) = (1-gamma(2))*F(2);
qin(4) = (1-gamma(1))*F(1);

% Levels
h = m./(rho*A);

% Outflows
qout = a.*sqrt(2*g*h);

% Differential equations (mass balances)
xdot = zeros(4,1);
xdot(1) = rho*(qin(1) + qout(3) - qout(1));
xdot(2) = rho*(qin(2) + qout(4) - qout(2));
xdot(3) = rho*(qin(3) - qout(3));
xdot(4) = rho*(qin(4) - qout(4));
end