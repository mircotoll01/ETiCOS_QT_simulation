c = parcluster('local');
c.NumWorkers = 6;
parpool(c, c.NumWorkers)
% 
% delete(gcp);

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
twell = 3.5*nm;                     % well thickness (m)
tbarr = 2.5*nm;                   % barrier thickness (m)
nbarr = 20;                        % number of barriers in the geometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NzPlot = 101;                     % number of plot points for the LDOS (if plotted)
NzSC   = 1001;                    % points for the staircase approximation (if used)
Vvet   = [0:0.005:0.3];            % vector of applied voltages (if used) (V)
Evet   = [-0.3501:1e-4:0.3]*eV;    % energy vector under study (J)
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

LDOS0 = zeros(length(Evet), NzPlot);

for indE = 1:length(Evet)
    E = Evet(indE);
    [GammaRight, Zinf, kz] = f_EvalGamma(E, Ueq, Leq, meffn);
    GammaLeft = f_EvalGamma(E, Ueq(end:-1:1), Leq(end:-1:1), meffn);
    GammaLeft = GammaLeft(end:-1:1);

    [LDOS0(indE,:), zvet] = f_EvalLDOS(NzPlot, GammaRight, GammaLeft, kz, Zinf, L);

    f0 = 1./(1+exp((E-0)./(kB*Temp)));

    eSpectralDensity(indE,:) = LDOS0(indE,:).*f0;
    eSpectralDensity(isnan(eSpectralDensity)) = 0;
end

eDensity0 = 1/(2*pi)*trapz(Evet, eSpectralDensity);

eDensity_rec_mean = mean(1./eDensity0(5:97));

for indV = 1:length(Vvet)
    T = zeros(length(Evet),1);
    SpectralJ = zeros(1, length(Evet));
    
    V = Vvet(indV);

    [U, L] = f_ApplyField(Ueq, Leq, V, NzSC);
    for indE = 1:length(Evet)
        E = Evet(indE);
        [GammaRight, Zinf, kz] = f_EvalGamma(E, U, L, meffn);
        GammaLeft = f_EvalGamma(E, U(end:-1:1), L(end:-1:1), meffn);
        GammaLeft = GammaLeft(end:-1:1);        
        tau = f_EvalTau(GammaRight, kz, L);
        T(indE) = abs(tau).^2;

        muL = -q*V; muR = 0;
        fL = 1./(1+exp((E-muL)./(kB*Temp)));
        fR = 1./(1+exp((E-muR)./(kB*Temp)));
        SpectralJ(indE) = -q/hbar*(fL-fR)*T(indE);
        SpectralJ(isnan(SpectralJ)) = 0;
        
    end

    J(indV) = 1/(2*pi)*trapz(Evet, SpectralJ); % from the Landauer-Büttiker formula
end 

% Converting L, U vectors in plottable vectors
[Uplot,Lplot] = f_Geom2Plot(U,L);

Efield = Vvet./Lplot(end);
eVelocity = J./q.*eDensity_rec_mean;

figure,
grid on
hold on
plot(Efield/(1/cm), eVelocity/cm, 'LineWidth',2)
xlabel('Electric field (V/cm)')
ylabel('Carrier velocity (cm/s)')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)

% this is the average electron velocity over the electric field
% Landau-Büttiker formula required for the exam
% this curve starts from zero because at zero field the current is zero since the difference of the Fermi functions (Fermi window) is zero (see Landau-Büttiker formula)
% current is zero at high fields because the transmission coefficient goes to zero. By applying an electric field the slope of the potential increases, and the transmission coefficient decreases. 
delete(gcp);