
function [Uplot,Lplot] = f_Geom2Plot(U,L)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Uplot,Lplot] = f_Geom2Plot(U,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input vectors U, L of the code are not ready to be plotted. The scope
% of this function is to convert them in quantity suitable to be pictured
% with the plot() command.
%
%
% Input parameters:
% o U     : vector with dimension equal to the number of layers, and
%           containing the value of the potential energy in each of them
%           (J)
% o L     : vector with dimension equal to the number of layers, and
%           containing the length of the layer in each of them (m)
%
% Output quantities:
% o Uplot : output vector of potential energies (J)
% o Lplot : output vector of lengths (m)
%
% In order to plot the geometry, it is possible to use the command, e.g.,
% plot(Lplot,Uplot)
% 
% Alberto Tibaldi, May 02, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uplot = zeros(1,2*length(U)-1);
Uplot(1:2:2*length(U)) = U;
Uplot(2:2:2*length(U)) = U;

Lplot = sort([0,cumsum(L),cumsum(L)+eps]);
Lplot = Lplot(1:end-1);
