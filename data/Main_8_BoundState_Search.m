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
nbarr = 1;                        % number of barriers in the geometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solver parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NzPlot = 101;                     % number of plot points for the LDOS (if plotted)
NzSC   = 1001;                    % points for the staircase approximation (if used)
Vvet   =  0:0.05:0.3;            % vector of applied voltages (if used) (V)
Evet   = (-0.1001:1e-3:0.3)*eV;    % energy vector under study (J)
% Evet   = 0.14*eV;                % energy vector under study (J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the geometry vectors U and L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = []; % initialization of the vector of the layer potentials
L = []; % initialization of the vector of the layer lengths
for ind = 1:nbarr
    U = [U, Ubarr,     0];
    L = [L, tbarr, twell];
end
% Terminate the geometry with a well
U = [U, Ubarr];
L = [L, tbarr];

for indE = 1:length(Evet)
    E = Evet(indE);
    zsrc = tbarr+0.5*twell;
    LoopGain(indE) = f_EvalLoopGain(E, zsrc, U, L, meffn);
end

[Uplot, Lplot] = f_Geom2Plot(U,L);

figure
set(gcf,'Position',[110 255 1154 420])
%
subplot(1,3,1)
grid on
hold on
box on
colormap(flipud(gray))
plot(Lplot/nm, Uplot/eV, 'r', 'LineWidth', 1.5)
axis([0,sum(L)/nm,min(Evet/eV),max(Evet/eV)])
xlabel('Position z (nm)')
ylabel('Potential energy (eV)')
set(gca,'FontName','Times New Roman','FontSize',12)
%
subplot(1,3,2)
grid on
hold on
box on
colormap(flipud(gray))
plot(abs(LoopGain), Evet/eV, 'r', 'LineWidth', 1.5)
axis([0,1.5,min(Evet/eV),max(Evet/eV)])
xlabel('|G_L|')
ylabel('Energy (eV)')
set(gca,'FontName','Times New Roman','FontSize',12)
%
subplot(1,3,3)
grid on
hold on
box on
colormap(flipud(gray))
plot(0*ones(size(Evet)), Evet/eV, 'k--', 'LineWidth', 1.5)
plot(angle(LoopGain)*180/pi, Evet/eV, 'r', 'LineWidth', 1.5)
axis([-200,200,min(Evet/eV),max(Evet/eV)])
xlabel('\angle G_L (deg)')
ylabel('Energy (eV)')
set(gca,'FontName','Times New Roman','FontSize',12)

fun = @(E)angle(f_EvalLoopGain(E, zsrc, U, L, meffn));
OPT = optimset;
OPT.TolX = 1e-10*eV;
E_BS = fzero(fun, 0.1*eV, OPT);

figure(1)
%eV
subplot(1,3,1)
II = ones(size(Lplot));
II([1,6]) = NaN;
plot(Lplot/nm,E_BS/eV*II,'m-.','LineWidth',1.5)
%
subplot(1,3,2)
plot(1,E_BS/eV,'mo','MarkerSize',8,'LineWidth',2)
%
subplot(1,3,3)
plot(0,E_BS/eV,'mo','MarkerSize',8,'LineWidth',2)