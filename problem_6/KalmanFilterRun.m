function out = KalmanFilterRun(kf,Y,U,xhat0,P0)
% Runs the static and the dynamic Kalman filter over a recorded data set
% (open loop: the inputs U do not depend on the estimates).
%
% In:  kf      design from KalmanFilterDesign
%      Y, U    measurements (ny x N) and inputs (nu x N), deviation variables
%      xhat0   initial one-step prediction x_{0|-1}
%      P0      initial covariance P_{0|-1} (dynamic filter)
% Out: out.XfS, out.XfD    filtered estimates x_{k|k}, static and dynamic
%      out.PfD             P_{k|k} of the dynamic filter (nx x nx x N)
%      out.nisS, out.nisD  normalized innovation squared e' Re^{-1} e
%      out.Kerr            ||K_k - Kfx|| / ||Kfx||

N  = size(Y,2);
nx = size(kf.A,1);
C  = kf.C;

out.XfS  = zeros(nx,N);
out.XfD  = zeros(nx,N);
out.PfD  = zeros(nx,nx,N);
out.nisS = zeros(1,N);
out.nisD = zeros(1,N);
out.Kerr = zeros(1,N);

xpS = xhat0;
xpD = xhat0;  PpD = P0;
for k = 1:N
    Kk          = (PpD*C')/(C*PpD*C' + kf.R);
    out.Kerr(k) = norm(Kk - kf.Kfx)/norm(kf.Kfx);

    [out.XfD(:,k),PfD,xpD,PpD,eD,ReD] = KalmanFilterDynamic(kf,xpD,PpD,Y(:,k),U(:,k));
    [out.XfS(:,k),xpS,eS]             = KalmanFilterStatic(kf,xpS,Y(:,k),U(:,k));

    out.PfD(:,:,k) = PfD;
    out.nisD(k)    = eD'*(ReD\eD);
    out.nisS(k)    = eS'*(kf.Re\eS);
end
end
