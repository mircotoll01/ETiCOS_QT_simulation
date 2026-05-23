function [tau] = f_EvalTau(GammaRight, kz, L)
    tau = 1;

    for ind = length(GammaRight)-1:-1:1
        tau = tau*exp(1i*kz(ind+1)*L(ind+1));
        GammaR = GammaRight(ind+1)*exp(2i*kz(ind+1)*L(ind+1));
        GammaL = GammaRight(ind);

        tau = tau*(1 + GammaL)/(1 + GammaR);
    end
end