function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)

% QPSOLVER  min 1/2 x'Hx + g'x  s.t.  l <= x <= u,  bl <= Ax <= bu 

% bl <= Ax <= bu  <=>  [A; -A] x <= [bu; -bl] in new form 

opts = optimoptions('quadprog', 'Display', 'off');

[x, fval, exitflag, output, lambda] = quadprog(H, g, [A; -A], [bu; -bl], [], [], l, u, xinit, opts);

info.fval     = fval;
info.exitflag = exitflag;
info.output   = output;
info.lambda   = lambda;

end