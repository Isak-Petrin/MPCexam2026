function y = QuadrupleTankSensor(x,p)
% Sensor output = measured tank levels [cm]

A   = p(5:8,1);
rho = p(12,1);

y = x./(rho*A);
end
