function kf = KalmanFilterDesign(A,B,C,G,Q,R)
% Design of a discrete-time (static/stationary) Kalman filter.
%
% System:
%   x_{k+1} = A*x_k + B*u_k + G*w_k
%   y_k     = C*x_k + v_k
%
%   w_k ~ N(0,Q),  v_k ~ N(0,R),  w and v independent (S = 0)
%
% The stationary one-step prediction covariance P = P_{k|k-1} solves the
% discrete-time algebraic Riccati equation (DARE)
%
%   P = A P A' + G Q G' - A P C' (C P C' + R)^{-1} C P A'
%
% It is computed without the Control System Toolbox (no dlqe/idare) by the
% structure-preserving doubling algorithm (Chu, Fan & Lin, 2005), which
% converges quadratically, i.e. in a few tens of iterations. Requires
% (A, C) detectable and (A, G Q^{1/2}) stabilizable.

nx = size(A,1);
Qw = G*Q*G';
Qw = (Qw + Qw')/2;

% --- Doubling algorithm for  P = A P (I + C'R^{-1}C P)^{-1} A' + Qw
Ak = A';
Gk = C'*(R\C);
Hk = Qw;
I  = eye(nx);
converged = false;
for it = 1:100
    W    = I + Gk*Hk;
    WiA  = W\Ak;                       % W^{-1} A_k
    WiG  = W\Gk;                       % W^{-1} G_k
    Hn   = Hk + Ak'*Hk*WiA;
    Gk   = Gk + Ak*WiG*Ak';
    Ak   = Ak*WiA;
    Gk   = (Gk + Gk')/2;
    Hn   = (Hn + Hn')/2;
    if norm(Hn - Hk,1) <= 1e-13*max(1,norm(Hn,1))
        Hk = Hn;
        converged = true;
        break
    end
    Hk = Hn;
end
P = Hk;

% Check the Riccati residual
Re  = C*P*C' + R;
res = A*P*A' + Qw - (A*P*C')/Re*(C*P*A') - P;
if ~converged || norm(res,1) > 1e-8*max(1,norm(P,1))
    warning('KalmanFilterDesign:DARE', ...
        'Riccati solution may be inaccurate (relative residual %.1e). Check detectability of (A,C).', ...
        norm(res,1)/max(1,norm(P,1)));
end

% Stationary measurement update
Kfx = P*C'/Re;

% Stationary filtered covariance P_{k|k}
Pf = P - Kfx*Re*Kfx';
Pf = (Pf + Pf')/2;

% Store model
kf.A  = A;
kf.B  = B;
kf.C  = C;
kf.Qw = Qw;
kf.R  = R;

kf.P   = P;
kf.Pf  = Pf;
kf.Re  = Re;
kf.Kfx = Kfx;

end
