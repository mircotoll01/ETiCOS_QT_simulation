
function [USC,LSC] = f_ApplyField(U,L,V,NzSC)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [USC,LSC] = f_ApplyField(U,L,V,NzSC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting from U, L in the usual format (e.g., that of f_EvalGamma), this
% function discretizes the profile in NzSC points (Number of z points for
% the StairCase). This allows to apply a constant electric field to the
% potential, corresponding to an additional linear potential starting from
% -qV applied to the left contact, ending to 0 at the right contact.
%
% Input parameters:
%
% o U     : vector with dimension equal to the number of layers, and
%           containing the value of the potential energy in each of them
%           (J)
% o L     : vector with dimension equal to the number of layers, and
%           containing the length of the layer in each of them (m)
% o V     : scalar, voltage to be applied at the left contact (V)
% o NzSC  : number of points for the staircase discretization
%
% Output quantities:
%
% o USC   : vector U after the staircase discretization and including the
%           constant electric field (J)
% o LSC   : vector L after the staircase discretization (m)
% 
% Alberto Tibaldi, May 02, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 1.6021766208e-019;  % elementary charge, C

% Converting U, L vectors in plottable vectors
[Uplot,Lplot] = f_Geom2Plot(U,L);

% Interpolate (linearly) the potential U, leading to UInterp.
zInterp = linspace(0,max(Lplot),NzSC);
UInterp = interp1(Lplot,Uplot,zInterp);

% Evaluating the potential energy associated to the electric field
Ufield = -q*linspace(V,0,NzSC);

% Summing UInterp and Ufield and performing staircase interpolation. Then,
% we add Ufield to UInterp.
[zp,Up] = stairs(zInterp,(UInterp+Ufield));

% Re-converting in vectors (Lnew,Unew) acceptable for the transmission line
% solver
zzp = zp(1:2:end);
UUp = Up(1:2:end);
LSC = diff(zzp).';
USC = UUp(1:end-1).';