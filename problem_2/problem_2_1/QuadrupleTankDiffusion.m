function G = QuadrupleTankDiffusion(~,~,~,~,psde)
% Diffusion matrix for the four-tank SDE model.
% psde = [p; sigma], where sigma = [sigma_F3; sigma_F4]

    rho     = psde(end-2);  % Last parameter in p
    sigma1  = psde(end-1);
    sigma2  = psde(end);

    G = zeros(4,2);
    G(3,1) = rho*sigma1;
    G(4,2) = rho*sigma2;
end