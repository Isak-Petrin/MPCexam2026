function m = Problem4StateSpace(Kid,Tid,gamma1,gamma2,Ts)
% Six-state MIMO realization of the transfer function model from Problem 4.4
%
%   G11 = K11/(T1 s+1)              G12 = K12/((T1 s+1)(T3 s+1))
%   G21 = K21/((T2 s+1)(T4 s+1))    G22 = K22/(T2 s+1)
%
% Each G_ij is realized in observer canonical form (Problem 4.6) and stacked:
%   x = [x1; x2; x3; x4],  x1, x4 in R,  x2, x3 in R^2
%   A = blkdiag(A1,A2,A3,A4),  B = [B1 0; 0 B2; B3 0; 0 B4],
%   C = [C1 C2 0 0; 0 0 C3 C4]
% Disturbances: F3 enters tank 3 as (1-gamma2)*F2 does and F4 enters tank 4
% as (1-gamma1)*F1 does, so F3 -> h1 = G12/(1-gamma2), F4 -> h2 = G21/(1-gamma1):
%   E = [0 0; B2/(1-gamma2) 0; 0 B3/(1-gamma1); 0 0]
%
% Returns the continuous-time (Ac,Bc,Cc,Ec) and ZOH discrete-time
% (Ad,Bd,Ed,C) matrices in deviation variables (u = F12-us, d = F34-ds, y = h12-hs).

% name, gain, denominator, input (1 = F1, 2 = F2), output, disturbance (0 = none)
blk = {'G11', Kid(1,1), [Tid(1) 1],                  1, 1, 0, 0;
       'G12', Kid(1,2), conv([Tid(1) 1],[Tid(3) 1]), 2, 1, 1, 1/(1-gamma2);
       'G21', Kid(2,1), conv([Tid(2) 1],[Tid(4) 1]), 1, 2, 2, 1/(1-gamma1);
       'G22', Kid(2,2), [Tid(2) 1],                  2, 2, 0, 0};

nb  = size(blk,1);
nu  = 2;  ny = 2;  nd = 2;
nxi = zeros(1,nb);
for i = 1:nb, nxi(i) = length(blk{i,3}) - 1; end
nx  = sum(nxi);
idx = [0 cumsum(nxi)];

Ac = zeros(nx);  Bc = zeros(nx,nu);  Cc = zeros(ny,nx);  Ec = zeros(nx,nd);
for i = 1:nb
    den = blk{i,3}/blk{i,3}(1);                % monic denominator
    b0  = blk{i,2}/blk{i,3}(1);                % constant numerator
    n   = nxi(i);
    a   = den(2:end);
    Ai  = [-a(:), [eye(n-1); zeros(1,n-1)]];   % observer canonical form
    Bi  = [zeros(n-1,1); b0];
    Ci  = [1 zeros(1,n-1)];
    r   = idx(i)+1:idx(i+1);
    Ac(r,r)        = Ai;
    Bc(r,blk{i,4}) = Bi;
    Cc(blk{i,5},r) = Ci;
    if blk{i,6} > 0
        Ec(r,blk{i,6}) = blk{i,7}*Bi;
    end
end

Phi = expm([Ac Bc Ec; zeros(nu+nd, nx+nu+nd)]*Ts);

m.Ac = Ac;  m.Bc = Bc;  m.Cc = Cc;  m.Ec = Ec;
m.Ad = Phi(1:nx, 1:nx);
m.Bd = Phi(1:nx, nx+1:nx+nu);
m.Ed = Phi(1:nx, nx+nu+1:end);
m.C  = Cc;
m.Ts = Ts;
end
