function [Ad,Qd] = DiscretizeSDE(A,G,Ts)
% Discretization of the process noise in  dX = A X dt + G dw  (Van Loan).
%
%   w_k ~ N(0,Qd),   Qd = int_0^Ts expm(A s) G G' expm(A' s) ds
%
%   expm([-A G*G'; 0 A'] Ts) = [Phi11 Phi12; 0 Phi22]
%   Ad = Phi22',  Qd = Phi22' * Phi12

nx = size(A,1);

L     = expm([-A, G*G'; zeros(nx), A']*Ts);
Phi12 = L(1:nx, nx+1:2*nx);
Phi22 = L(nx+1:2*nx, nx+1:2*nx);

Ad = Phi22';
Qd = Ad*Phi12;
Qd = (Qd + Qd')/2;   % remove round-off asymmetry
end
