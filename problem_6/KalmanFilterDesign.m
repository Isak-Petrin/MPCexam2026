function kf = KalmanFilterDesign(A,B,C,G,Q,R)
% Design of the Kalman filter for the discrete-time linear system
%   x_{k+1} = A x_k + B u_k + G w_k,   w_k ~ N(0,Q)
%   y_k     = C x_k + v_k,             v_k ~ N(0,R)
%
% Stores the model for the dynamic filter (KalmanFilterDynamic) and
% computes the stationary gain for the static filter (KalmanFilterStatic)
% with dlqe (Control System Toolbox), which solves the discrete algebraic
% Riccati equation
%   P = A P A' + G Q G' - A P C' (C P C' + R)^-1 C P A'
%
% For offset-free estimation, pass the disturbance-augmented model.

[Kfx,P,Pf] = dlqe(A,G,C,Q,R);

kf.A = A; kf.B = B; kf.C = C; kf.Qw = G*Q*G'; kf.R = R;
kf.P   = P;              % stationary P_{k|k-1}
kf.Re  = C*P*C' + R;     % stationary innovation covariance
kf.Kfx = Kfx;            % stationary filter gain, P*C'/Re
kf.Pf  = Pf;             % stationary P_{k|k}
end
