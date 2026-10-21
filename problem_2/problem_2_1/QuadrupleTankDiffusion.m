function G = QuadrupleTankDiffusion(p, sigma1, sigma2)

% træk rho ud af parameter vektor p

    G = zeros(4,2);
    G(3,1) = rho * sigma1;
    G(4,2) = rho * sigma2;
end
