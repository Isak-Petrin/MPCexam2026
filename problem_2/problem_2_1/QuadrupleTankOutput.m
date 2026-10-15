function z = QuadrupleTankOutput(x,p)
% Controlled outputs z = h(x,p): lower tank levels [cm]

A   = p(5:8,1);
rho = p(12,1);

h = x./(rho*A);
z = h(1:2);
end
