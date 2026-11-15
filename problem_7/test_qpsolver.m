%% test_qpsolver.m - Problem 7: tests of the QP solver interface
%
%   [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)
%
%   min  1/2 x'Hx + g'x   s.t.  l <= x <= u,  bl <= Ax <= bu
%
% Every case has a solution known in advance, so the tests do not depend on
% the multipliers returned by quadprog. Cases 1-7 use the 2-variable QP
%
%   min (x1-1)^2 + (x2-2.5)^2      (H = 2I, g = [-2;-5])
%
% whose unconstrained minimiser is (1, 2.5). Case 8 is a random strictly
% convex QP built backwards from a chosen solution and chosen multipliers
% (KKT conditions), so its solution is unique and known exactly.
%
% Requires the Optimization Toolbox (quadprog).

clear; clc;

tol = 1e-5;                       % tolerance on ||x - xref||_inf
H2  = 2*eye(2);
g2  = [-2; -5];
x0  = zeros(2,1);
big = 10*ones(2,1);

cases = {};

% 1) No active constraints -> unconstrained minimiser -H\g
cases{end+1} = struct('name','1 inactive constraints', ...
    'H',H2,'g',g2,'l',-big,'u',big,'A',[1 1],'bl',-100,'bu',100, ...
    'x0',x0,'xref',[1; 2.5],'feasible',true);

% 2) Bounds only, empty A, upper bound on x2 active
cases{end+1} = struct('name','2 active bound, A empty', ...
    'H',H2,'g',g2,'l',-big,'u',[10; 2],'A',zeros(0,2),'bl',zeros(0,1),'bu',zeros(0,1), ...
    'x0',x0,'xref',[1; 2],'feasible',true);

% 3) General constraint active at its upper side: x1 + x2 <= 2
cases{end+1} = struct('name','3 active bu', ...
    'H',H2,'g',g2,'l',-big,'u',big,'A',[1 1],'bl',-10,'bu',2, ...
    'x0',x0,'xref',[0.25; 1.75],'feasible',true);

% 4) General constraint active at its lower side: x1 + x2 >= 5
cases{end+1} = struct('name','4 active bl', ...
    'H',H2,'g',g2,'l',-big,'u',big,'A',[1 1],'bl',5,'bu',10, ...
    'x0',x0,'xref',[1.75; 3.25],'feasible',true);

% 5) Equality constraint written as bl = bu: x1 + x2 = 2
cases{end+1} = struct('name','5 equality bl = bu', ...
    'H',H2,'g',g2,'l',-big,'u',big,'A',[1 1],'bl',2,'bu',2, ...
    'x0',x0,'xref',[0.25; 1.75],'feasible',true);

% 6) One-sided constraint and free variables: infinite l, u and bl
cases{end+1} = struct('name','6 infinite l, u, bl', ...
    'H',H2,'g',g2,'l',-inf(2,1),'u',inf(2,1),'A',[1 1],'bl',-inf,'bu',2, ...
    'x0',x0,'xref',[0.25; 1.75],'feasible',true);

% 7) Infeasible: 0 <= x <= 1 but x1 + x2 >= 5
cases{end+1} = struct('name','7 infeasible', ...
    'H',H2,'g',g2,'l',[0; 0],'u',[1; 1],'A',[1 1],'bl',5,'bu',10, ...
    'x0',x0,'xref',[NaN; NaN],'feasible',false);

% 8) Random strictly convex QP (n = 20, m = 8) with a constructed solution.
%    Choose x*, the active set and positive multipliers, then set g from
%    stationarity:  H x* + g - lam_l + lam_u + A'(mu_u - mu_l) = 0
rng(1);
n = 20; m = 8;
M  = randn(n);
H  = M'*M + eye(n);
A  = randn(m,n);
xs = randn(n,1);
l  = xs - 1;   u  = xs + 1;
bl = A*xs - 1; bu = A*xs + 1;
lam_l = zeros(n,1); lam_u = zeros(n,1);
mu_l  = zeros(m,1); mu_u  = zeros(m,1);
l(1:3)  = xs(1:3);    lam_l(1:3) = [1; 2; 3];     % active lower bounds
u(4:6)  = xs(4:6);    lam_u(4:6) = [1; 2; 3];     % active upper bounds
bu(1:2) = A(1:2,:)*xs; mu_u(1:2) = [1; 2];        % active at bu
bl(3:4) = A(3:4,:)*xs; mu_l(3:4) = [1; 2];        % active at bl
g = -H*xs + lam_l - lam_u - A'*(mu_u - mu_l);
cases{end+1} = struct('name','8 random QP n=20 m=8', ...
    'H',H,'g',g,'l',l,'u',u,'A',A,'bl',bl,'bu',bu, ...
    'x0',zeros(n,1),'xref',xs,'feasible',true);

%% Run
nc = numel(cases);
pass = false(nc,1);
fprintf('%-26s %-5s %9s %11s %11s  %s\n', ...
    'case','res','exitflag','err_x','viol','note');
fprintf('%s\n', repmat('-',1,84));

for i = 1:nc
    c = cases{i};
    note = '';
    err = NaN; viol = NaN; ef = NaN;
    lastwarn('');
    try
        [x,info] = qpsolver(c.H,c.g,c.l,c.u,c.A,c.bl,c.bu,c.x0);
        ef = info.exitflag;
        if c.feasible
            err  = norm(x - c.xref, inf);
            Ax   = c.A*x;
            viol = max([0; c.l - x; x - c.u; c.bl - Ax; Ax - c.bu]);
            pass(i) = ef > 0 && err <= tol && viol <= tol;
        else
            pass(i) = ef < 0;          % solver must report failure
        end
        w = lastwarn;
        if ~isempty(w), note = ['warning: ' w]; end
    catch ME
        note = ['ERROR: ' ME.message];
    end
    res = 'FAIL'; if pass(i), res = 'PASS'; end
    fprintf('%-26s %-5s %9g %11.2e %11.2e  %s\n', c.name, res, ef, err, viol, note);
end

fprintf('%s\n', repmat('-',1,84));
fprintf('%d of %d cases passed (tol = %.0e)\n', sum(pass), nc, tol);
