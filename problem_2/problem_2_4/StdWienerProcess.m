function [W,Tw,dW] = StdWienerProcess(T,N,nW,Ns,seed)
% Multivariate standard Wiener process  (Lecture 03B slide 6)
% W  : nW x (N+1) x Ns  realizations, W(:,1,:) = 0
% Tw : 1 x (N+1)        time grid
% dW : nW x N x Ns      increments

    if nargin == 5
        rng(seed);
    end

    dt = T/N;

    dW = sqrt(dt)*randn(nW,N,Ns);
    W  = [zeros(nW,1,Ns) cumsum(dW,2)];

    Tw = 0:dt:T;
end
