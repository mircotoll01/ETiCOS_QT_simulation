function LoopGain = f_EvalLoopGain(E, zsrc, U, L, meffn)
    [GammaRight, Zinf, kz] = f_EvalGamma(E, U, L, meffn);

    GammaLeft = f_EvalGamma(E, U(end:-1:1), L(end:-1:1), meffn);
    GammaLeft = GammaLeft(end:-1:1);

    zint = cumsum([0 L]);

    indLayer = find(zsrc >= zint);
    indLayer = indLayer(end);

    zint_Right = sum([L(1:indLayer)]);
    GammaRightSRC = GammaRight(indLayer)*exp(2i*kz(indLayer)*abs(zsrc-zint_Right));

    zint_Left = sum([L(1:indLayer-1)]);
    GammaLeftSRC = GammaLeft(indLayer)*exp(2i*kz(indLayer)*abs(zsrc-zint_Left));

    LoopGain = GammaRightSRC*GammaLeftSRC;
end