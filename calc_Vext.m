function [Vext] = calc_Vext(p_electrode, nodes_rnv, Istim)
%Vext 50mV ~for activation

rho_e = 350 / 100; %350 ohm * cm * 1 m / 100 cm \cite{major2005}

n_distance = pdist2(p_electrode, nodes_rnv);

Vext = rho_e .* Istim ./ (4 .* pi .* n_distance);
