clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unit conversion constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The engine of the code works with SI units; in this view, conversion
% constants are useful to insert numbers in more "confortable", yet less
% standard, units
cm = 1e-2;                        % multiply to convert from cm to m
um = 1e-6;                        % multiply to convert from um to m
nm = 1e-9;                        % multiply to convert from nm to m
eV = 1.6021766208e-019;           % multiply to convert from eV to J

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q    = 1.6021766208e-019;         % elementary charge (C)
m0   = 9.10938291e-31;            % electron mass (kg)
kB   = 1.3806488e-23;             % Boltzmann constant (J/K)
h    = 6.626070040e-34;           % Planck constant (J*s)
hbar = h/(2*pi);                  % reduced Planck constant (J*s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Temp  = 300;                      % temperature (K)
meffn = 0.067*m0;                 % effective mass (kg)
Ubarr = 0.25*eV;                  % barrier potential (J)
twell = 3*nm;                     % well thickness (m)
tbarr = 2.5*nm;                   % barrier thickness (m)
nbarr = 20;                        % number of barriers in the geometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NzPlot = 101;                     % number of plot points for the LDOS (if plotted)
NzSC   = 1001;                    % points for the staircase approximation (if used)
Vvet   = [0:0.05:0.3];            % vector of applied voltages (if used) (V)
Evet   = [-0.1001:2e-4:0.25]*eV;    % energy vector under study (J)
%Evet   = 0.14*eV;                % energy vector under study (J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the geometry vectors U and L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = []; % initialization of the vector of the layer potentials
L = []; % initialization of the vector of the layer lengths
for ind = 1:nbarr
    U = [U, 0,     Ubarr];
    L = [L, twell, tbarr];
end
% Terminate the geometry with a well
U = [U,     0];
L = [L, twell];

Ueq = U;
Leq = L;

V = 0.1;

[U, L] = f_ApplyField(Ueq, Leq, V, NzSC);
for indE = 1:length(Evet)
    E = Evet(indE);
    [GammaRight, Zinf, kz] = f_EvalGamma(E, U, L, meffn);
    GammaLeft = f_EvalGamma(E, U(end:-1:1), L(end:-1:1), meffn);
    GammaLeft = GammaLeft(end:-1:1);

    [LDOS(indE,:), zvet] = f_EvalLDOS(NzPlot, GammaRight, GammaLeft, kz, Zinf, L);
    
    tau = f_EvalTau(GammaRight, kz, L);
    T(indE) = abs(tau).^2;

    muL = -q*V; muR = 0;
    fL = 1./(1+exp((E-muL)./(kB*Temp)));
    fR = 1./(1+exp((E-muR)./(kB*Temp)));
    SpectralJ(indE) = -q/hbar*(fL-fR)*T(indE);
    SpectralJ(isnan(SpectralJ)) = 0;
    
end

J = 1/(2*pi)*trapz(Evet, SpectralJ)

% Converting L, U vectors in plottable vectors
[Uplot,Lplot] = f_Geom2Plot(U,L);
% Plotting local density of states
figure
%
set(gcf,'Position',[110 255 1154 420])
subplot(1,2,1)
grid on
hold on
box on
plot(T,Evet/eV,'b','LineWidth',1.5)
ylim([min(Evet/eV),max(Evet/eV)])
xlabel('Transmittance')
ylabel('Energy (eV)')
title('Transmittance spectrum')
set(gca,'FontName','Times New Roman','FontSize',12)
%
subplot(1,2,2)
grid on
hold on
box on
colormap(flipud(gray))
imagesc(zvet/nm,Evet/eV, real(LDOS/(1/eV)/(1/cm)))
plot(Lplot/nm, Uplot/eV, 'b', 'LineWidth', 1.5)
axis([min(zvet/nm),max(zvet/nm),min(Evet/eV),max(Evet/eV)])
xlabel('Position z (nm)')
ylabel('Energy (eV)')
title('LDOS (1/(eV cm))')
colorbar
set(gca,'FontName','Times New Roman','FontSize',12)
