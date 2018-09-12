function [activity_bool] = is_activated(mu_nodes,neuropoints, p_elec,t_pulse, Istim)

if nargin == 3
    t_pulse = 0.5e-3; %pulsewidth in (s)
    Istim = 10e-3; %stimulus in (A)
end

npnew2 = subset_tree(neuropoints, mu_nodes);
[mu_pts,dtrs] = dot_tree(npnew2);
vxt = calc_Vext(p_elec/100, mu_pts/100, Istim);
adj_mat = create_adj_matrix(dtrs, npnew2);

D = 10e-6;

activity = neuron_response(vxt',adj_mat,D,t_pulse);
activity_bool = any(activity);