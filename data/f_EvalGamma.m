function [GammaRight, Zinf, kz] = f_EvalGamma(E, U, L, meffn)
    hbar = 6.626070040e-34/(2*pi);

    kz = sqrt(2*meffn./(hbar^2)*(E-U));

    Yinf = sqrt(8./meffn.*(E-U));
    Zinf = 1./Yinf;

    GammaRight = zeros(size(U));

    for ind = length(GammaRight)-1:-1:1
        GammaRp = GammaRight(ind+1)*exp(2i*kz(ind+1)*L(ind+1));
        GammaLR = (Zinf(ind+1) - Zinf(ind))/(Zinf(ind+1) + Zinf(ind));
        GammaRight(ind) = (GammaLR + GammaRp)/(1 + GammaLR*GammaRp);
    end

end