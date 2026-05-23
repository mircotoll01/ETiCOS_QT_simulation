function [LDOS, zvet] = f_EvalLDOS(NzPlot, GammaRight, GammaLeft, kz, Zinf, L)
    hbar = 6.626070040e-34/(2*pi);
    zint = cumsum([0 L]);
    iring = -1i*4/hbar;

    zvet = linspace(0, max(zint)*(1-1e-10), NzPlot);

    for indz = 1:NzPlot
        zsrc = zvet(indz);
        indLayer = find(zsrc >= zint);
        indLayer = indLayer(end);

        zint_Right = sum([L(1:indLayer)]);
        GammaRightSRC = GammaRight(indLayer)*exp(2i*kz(indLayer)*abs(zsrc-zint_Right));
        ZRightSRC = Zinf(indLayer)*(1 + GammaRightSRC)/(1 - GammaRightSRC);

        zint_Left = sum([L(1:indLayer-1)]);
        GammaLeftSRC = GammaLeft(indLayer)*exp(2i*kz(indLayer)*abs(zsrc-zint_Left));
        ZLeftSRC = Zinf(indLayer)*(1 + GammaLeftSRC)/(1 - GammaLeftSRC);

        V(indz) = ZLeftSRC*ZRightSRC/(ZLeftSRC + ZRightSRC) * iring;

    end 

    LDOS = 1i*(V - conj(V));

end