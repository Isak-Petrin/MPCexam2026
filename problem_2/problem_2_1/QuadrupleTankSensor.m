function y = QuadrupleTankSensor(x,p)
% Sensor output y = g(x,p): measured levels in tank 1 and tank 2 [cm]

A   = p(5:8,1);
rho = p(12,1);

h = x./(rho*A);
y = h(1:2);
end
