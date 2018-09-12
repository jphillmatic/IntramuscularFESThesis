nMu = 40;
nNMJ = nMu *5;
LoadColorScales
global scale_discrete
global scale_gradient
global bloo
global reed

soleus_length = 39.3; %cm
soleus_volume = 356+58; %cm3
soleus_radius = sqrt(soleus_volume / soleus_length / pi); %cm

n = 5;
nminorbranches = 0;

unbalanced1 = true;
unbalanced2 = true;

while (nminorbranches < 18) || unbalanced1 || unbalanced2
    np0 = FraktalTwigTree872018(20);
    [x0, y0, r, ht] = find_cyl(np0);
    
    np0 = sortrows(np0,6);
    np0(:,3) = np0(:,3);
    ht_mod = soleus_length/ht;
    np0(:,3) = np0(:,3)*ht_mod;
    
    [tta rda] = cart2pol(np0(:,1), np0(:,2));
    rda = rda * soleus_radius / r;
    [np0(:,1), np0(:,2)] = pol2cart(tta,rda);
    
    np1 = removeOutliers(np0);
    
    nminorbranches = length(find_nmj(np1));
    n_nodes = max(np1(:,6));
    n_subbranches = find_nmj(np1);
    
    [main_branches, n_main_branches] = main_branch_nmjs(np1);
    
    balance_mat = [];
    for i = 1:size(n_main_branches)
        balance_mat(:,i) = n_main_branches(i)./n_main_branches;
    end
    if any(any(balance_mat < 0.5))
        unbalanced1 = true;
    else
        unbalanced1 = false;
    end
    
    [x0, y0, r, ht] = find_cyl(np1);

    rnmj = rand(nNMJ,1)*r*1.0;
    thetanmj = rand(nNMJ,1)*2*pi;
    znmj = rand(nNMJ,1)*3*ht/10 + ht/5 ;
    
    [xnmj, ynmj, znmj] = pol2cart(thetanmj,rnmj,znmj);
    nmj_pos = [xnmj(:)+x0 ynmj(:)+y0 znmj(:)];
    neuropoints = find_nearest_nodes(nmj_pos, np1);
    
    [main_branches, n_main_branches] = main_branch_nmjs(neuropoints);
    
    balance_mat = [];
    for i = 1:size(n_main_branches)
        balance_mat(:,i) = n_main_branches(i)./n_main_branches;
    end
    if any(any(balance_mat < 0.2))
        unbalanced2 = true;
    else
        unbalanced2 = false;
    end
end


[x0, y0, r, ht] = find_cyl(neuropoints);
figure(1)
clf
scatter3(nmj_pos(:,1),nmj_pos(:,2), nmj_pos(:,3), 'r')
plot_tree(neuropoints, [], bloo)
% plot_tree(neuropoints, [], bloo)
% scatter3(neuropoints(n_subbranches,1), neuropoints(n_subbranches,2), neuropoints(n_subbranches,3), 'r')
viscircles([x0,y0],r);
viscircles([x0,y0],r*0.8);
sl = ht/3;%r*1.2;
axis_vals = [-sl sl -sl sl 0 ht];
axis(axis_vals)
view(-55, 20)
view(-0, 90)

% %% Gif plot
% fname = 'spinMU';
% spin_plot(100, fname)
sizemod = 1.5;
set(gcf, 'Position', [0, 0, 600*sizemod, 600*sizemod])

% pause;
%% Tree subsetting
%Time to peak (P) and max twitch force(T) from Fugelvand 1993  
forceRange = 120; % 10-fold force diffecence between largest and smallest MU
fatigueRange = 180; % range of fatigue rates across the motor units (300 best)
contractionRange = 3; % Range in contraction times,
ct_max = 90;

[P, fatigue, arborLengths, m_csa]  = define_motor_pool(nMu, forceRange, fatigueRange, contractionRange);
tree_length = getTreeLength(neuropoints);

% n_arborLengths = ceil(arborLengths * 0.15 * n_nodes);
n_arborLengths = ceil(1.5*arborLengths / sum(arborLengths) * tree_length);

mu_tree = select_MU_nodes(n_arborLengths, neuropoints);
% mu_tree = simple_MU_nodes(nMu, neuropoints); 
npnew = subset_tree(neuropoints, unique([mu_tree{:}]));

[mu_thru] = findPathWeights(npnew, mu_tree);
mu_thru(:,2) = (mu_thru(:,2) / max(mu_thru(:,2))).^0.5;
clf
% plot_tree(neuropoints, [],scale_discrete(9,:));
for i = 1:nMu
    pmu = nMu + 1 - i ;
%     create_nRv(neuropoints, mu_tree{i});
    plot_tree(npnew, mu_tree{pmu},scale_gradient(i,:), mu_thru);
end
title(' ')
whitebg('black')
set(gcf, 'Color', [0 0 0])
% whitebg('white')
% set(gcf, 'Color', [0 0 0])
set(gcf, 'Position', [0, 0, 600*1.5, 600*1.5])
[x0, y0, r, ht] = find_cyl(npnew);
sl = ht/4;
axis_vals = [-sl sl -sl sl 0 ht];
axis(axis_vals)

view(-30, 25)
% %%
% % export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\SubTree.png' -transparent -m3
% % spin_plot(200,[-sl sl -sl sl 0 h],'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\spin_tree_white')
% % viscircles([x0,y0],r)
% 
% %% Treecutting
% % 
% % for mu = 1:nMu
% %     clf
% %     plot_tree(npnew, [], scale_discrete(9,:), mu_thru);
% %     plot_tree(npnew, mu_tree{mu}, scale_gradient(mu,:), mu_thru);
% %     axis([-sl sl -sl sl 0 h])
% %     axis vis3d
% %     view(-55, 15)
% %     pause;
% % end
% 
% % cycle_mu_plot(npnew,mu_tree,axis_vals,'\mu_cycle')
spin_plot(150, [-sl sl -sl sl 0 ht], 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Project Files\Jake Model\Figures\spinmus')
% % export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\Ding2002Force.png' -transparent -m3
%% Branch Property testing
for mu = 1:nMu
    npmu = subset_tree(neuropoints,mu_tree{mu});
    for i = 1:length(npmu(:,6)) 
        nd(i,mu) = length(find_daughters(npmu(i,6), npmu));
    end
end
    
tot_brnchs = sum(sum(nd>1));
brnchs_rate = sum(sum(nd>1))/nMu
tot_xb = [sum(sum(nd==2))/tot_brnchs, sum(sum(nd==3))/tot_brnchs, sum(sum(nd==4))/tot_brnchs, sum(sum(nd==5))/tot_brnchs]
%% Treecutting

for mu = 1:nMu
    clf
    plot_tree(npnew, [], scale_discrete(9,:), mu_thru);
    plot_tree(npnew, mu_tree{mu}, scale_gradient(mu,:), mu_thru);
    axis([-sl sl -sl sl 0 ht])
    axis vis3d
    view(-55, 15)
    pause;
end
%% NRV sequencing
% for i = 1:length(adj_mat)
%     clf
%     view(-45, 25)
%     axis([-sl sl -sl sl 0 h])
%     pnds = find(adj_mat(i,:));
%     plot_tree(npnew,mu_tree{mu}, scale_discrete(1,:));
%     scatter3(mu_pts(pnds,1),mu_pts(pnds,2),mu_pts(pnds,3))
%     scatter3(mu_pts(i,1),mu_pts(i,2),mu_pts(i,3),'r')
%     pause
% end

%% Electrode Distances
% ~15 minutes
%Points in motor unit 1
figure(9)
clf
hold on
% tree1pts = dot_tree(neuropoints);
% tree1pts = dot_tree(npnew);
% plot_tree(tree1pts, [], bloo)
% dot_tree(neuropoints,mu_tree{2})

%meshgrid setup
THETA = linspace(0, 2*pi, 9);
THETA = THETA(1:end-1);
RHO = linspace(r*0.1, r*1.25, 5);
Z = linspace(0,ht,5);

%Create meshgrid of points in polar coordinates
[theta,rho,z] = meshgrid(THETA,RHO,Z);

%flatten
o_pts = [theta(:), rho(:), z(:)];
%remove duplicates (rho = 0)
o_pts = unique(o_pts,'Rows');
o_pts = sortrows(o_pts,[3,2,1]);
%Redefine theta rho z
theta = o_pts(:,1); rho = o_pts(:,2); z = o_pts(:,3);


%convert polar coordinates to cartesian
[elec_x , elec_y, elec_z] = pol2cart(theta, rho, z);
%flatten cartesian coords
e_pts = [elec_x(:), elec_y(:), elec_z(:)];
% e_pts = [r,0,h/2];

%calculate distances from each cartesian point to each point in th tree

% dsts = pdist2(e_pts,tree1pts(:,1:3));
percent_force_recruitment = zeros(size(e_pts,1),1);
active_mu_bool = logical(zeros(size(e_pts,1),nMu));
np_mesh_activated = cell(nMu,1);
tic
for i = 1:size(e_pts,1)
    
    p_elec = e_pts(i,:);
    
    for mu = 1:nMu
        active_mu_bool(i, mu) = is_activated(mu_tree{mu}, npnew, p_elec, 0.5e-3, 10e-3 * 0.7);
    end
    percent_force_recruitment(i) = sum(P(active_mu_bool(i,:)))/sum(P);
    
    np_mesh_activated{i} = unique([mu_tree{find(active_mu_bool(i,:))}]);
%     if ~isempty(a)
%         np_activated{i} = a;
%     else
%         np_activated{i} = 0;
%     end
    %     [mu_pts,dtrs] = dot_tree(npnew);
    %     vxt = calc_Vext(p_elec/100, mu_pts/100, 10);
    %     adj_mat = create_adj_matrix(dtrs, npnew);
    %     D = 10e-6;
    %     t_pulse = 0.5e-3;
    %     activity = neuron_response(vxt'*1e-3,adj_mat,D,t_pulse);
    
%     np_activated = unique([mu_tree{find(active_bool)}]);
    
    %     figure(10)
    %     clf
    %     hold on
    %     grid on
    %
    %     plot_tree(npnew, [], scale_discrete(1,:));
    %     plot_tree(npnew, np_activated, scale_discrete(2,:));
    %
    %     %     pause(3)
    %     %     plot3(mu_pts(find(activity),1),mu_pts(find(activity),2),mu_pts(find(activity),3),'.','Color', reed)
    %     %     plot3(mu_pts(~activity,1),mu_pts(~activity,2),mu_pts(~activity,3),'.','Color', bloo)
    %     plot3(p_elec(1),p_elec(2),p_elec(3),'o','Color', scale_discrete(8,:))
    %     plot3(x0+r_pos,repmat(p_elec(2),1,steps),repmat(p_elec(3),1,steps),'o','Color', scale_discrete(9,:))
    %     axis([-sl sl -sl sl 0 h])
    %     view(10, 11)
    
end
%%
% [pt2treemin,ys] = min(dsts,[],2);

% scatter3(e_pts(:,1), e_pts(:,2), e_pts(:,3), '.', 'MarkerEdgeAlpha',0.75)
% plot_tree(neuropoints, [], bloo)
%get groups by rho and z
groupd = findgroups(rho,z);
%apply mean function to distances by groups of rho-z combinations
gmeans = splitapply(@mean,percent_force_recruitment,groupd);
distnc = flip(reshape(gmeans, length(Z), length(RHO)),1);

[rr, elec_z] = meshgrid(RHO,Z);
scatter3(rr(:),zeros(length(rr(:)),1), elec_z(:), flip(distnc(:),1)*100+1,'r')
hold off
view([0 0])
%%
cycle_activation_plot2(npnew,mu_tree,active_mu_bool,e_pts, axis_vals,percent_force_recruitment','CycleActivations2.0')
%%
cycle_position_plot(npnew,mu_tree,active_mu_bool,e_pts, axis_vals, percent_force_recruitment', 'CycleActivations2.0')
%%
i_stims = linspace(2.5,20, 6)*1e-3;
cycle_current_plot(npnew,mu_tree,i_stims, 'Cycle Current')
%%
% distnc(1) = 1;
figure(10)
clf
rmp = RHO./r;
zmp = flip(Z./ht);
h = heatmap(rmp,zmp, distnc,'CellLabelColor','k', 'CellLabelFormat','%0.2g');
h.Title = 'Mean Proportion of Total Force Activated';
h.XLabel = 'Normalized Radial Position';
h.YLabel = 'Normalized Proximal-Distal Position';
