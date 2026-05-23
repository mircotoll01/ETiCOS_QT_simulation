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
nbarr = 2;                        % number of barriers in the geometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NzPlot = 101;                     % number of plot points for the LDOS (if plotted)
NzSC   = 1001;                    % points for the staircase approximation (if used)
Vvet   = [0:0.05:0.3];            % vector of applied voltages (if used) (V)
%Evet   = [0.0001:1e-3:1.2]*eV;    % energy vector under study (J)
Evet   = 0.14*eV;                % energy vector under study (J)

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

[GammaRight, Zinf, kz] = f_EvalGamma(Evet, U, L, meffn);

GammaLeft = f_EvalGamma(Evet, U(end:-1:1), L(end:-1:1), meffn);
GammaLeft = GammaLeft(end:-1:1);

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

figure(1)
plot(zvet/nm, LDOS, 'LineWidth', 2)
xlabel('Position (nm)')
ylabel('LDOS')
title('Local Density of States')
grid on
