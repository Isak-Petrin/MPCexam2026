function [Ad,Bd,Ed] = DiscretizeZOH(A,B,E,Ts)
% Exact zero-order-hold discretization of  Xdot = A X + B U + E D.
%
%   X_{k+1} = Ad X_k + Bd U_k + Ed D_k
%
% with Ad = expm(A Ts) and [Bd Ed] = int_0^Ts expm(A s) ds [B E].
% Computed with one matrix exponential of the augmented system in which
% U and D are extra states with zero dynamics (ZOH):
%
%   expm([A B E; 0 0 0] Ts) = [Ad Bd Ed; 0 I 0; 0 0 I]

nx = size(A,1);
nu = size(B,2);
nd = size(E,2);

M   = [A B E; zeros(nu+nd, nx+nu+nd)];
Phi = expm(M*Ts);

Ad = Phi(1:nx, 1:nx);
Bd = Phi(1:nx, nx+1:nx+nu);
Ed = Phi(1:nx, nx+nu+1:nx+nu+nd);
end
