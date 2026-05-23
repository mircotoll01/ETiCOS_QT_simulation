
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
tbarr = 2.5*nm;                     % barrier thickness (m)
nbarr = 1;                        % number of barriers in the geometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NzPlot = 101;                     % number of plot points for the LDOS (if plotted)
NzSC   = 1001;                    % points for the staircase approximation (if used)
Vvet   = [0:0.05:0.3];            % vector of applied voltages (if used) (V)
%Evet   = [0.0001:1e-3:1.2]*eV;  % energy vector under study (J)
Evet   = 0.14*eV;                  % energy vector under study (J)

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

E = Evet;
kz = sqrt(2*meffn./(hbar^2)*(E-U));

Yinf = sqrt(8./meffn.*(E-U));
Zinf = 1./Yinf;

GammaRight = zeros(size(U));

for ind = length(GammaRight)-1:-1:1
    GammaRp = GammaRight(ind+1)*exp(2i*kz(ind+1)*L(ind+1));
    GammaLR = (Zinf(ind+1) - Zinf(ind))/(Zinf(ind+1) + Zinf(ind));
    GammaRight(ind) = (GammaLR + GammaRp)/(1 + GammaLR*GammaRp);
end

GammaRight(1)

tau = 1;

for ind = length(GammaRight)-1:-1:1
    tau = tau*exp(1i*kz(ind+1)*L(ind+1));
    GammaR = GammaRight(ind+1)*exp(2i*kz(ind+1)*L(ind+1));
    GammaL = GammaRight(ind);

    tau = tau*(1 + GammaL)/(1 + GammaR);
end

tau

