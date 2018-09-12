%% TODO: MAKE TREES IN 3D SPACES
%%
%%Initialize values
%nodeData %
%1. x0
%2. y0
%3. z0
%4. level
%5. branch #
%6. node
%7. Parent
global scale_discrete
global scale_gradient
nMu = 60;
n = 10;
dmu = nMu*2;
mu=0;
inc = 0;

%Generate a fractal tree with at least x top level nodes
while mu < dmu
    neuropoints = FraktalTree(n,1,1,1);
    nmj_nodes = find_nmj(neuropoints);
    mu = length(nmj_nodes);
%     n_mu = length(find(neuropoints(:,4) == n));
    inc = inc + 1;
    if inc > n
        n = n + 1;
        inc = 0;
    end
end

% %generate fractal tree with points in a cylinder
% while n_mu < 120
%     neuropoints = FraktalTree(n,1,1,1);
%     n_mu = length(find(sum(neuropoints(:,1:2).^2,2) < 1));    inc = inc + 1;
%     if inc > 3
%         n = n + 1;
%         inc = 0;
%     end
% end

% neuropoints = neuropoints((sum(neuropoints(:,1:2).^2,2) < 1),:);
neuropoints = sortrows(neuropoints,6);
neuropoints(:,3) = neuropoints(:,3)*.75;
n_nodes = max(neuropoints(:,6));
figure(1)
clf
plot_tree(neuropoints, [], scale_discrete(1,:))
scatter3(neuropoints(nmj_nodes,1), neuropoints(nmj_nodes,2), neuropoints(nmj_nodes,3), 'r')
[x0, y0, r, h] = find_cyl(neuropoints);

sl = h/2;%r*1.2;
axis_vals = [-sl sl -sl sl 0 h];
axis(axis_vals)
view(-55, 20)
% %% Gif plot
% fname = 'spinMU';
% spin_plot(100, fname)
sizemod = 1.5;
set(gcf, 'Position', [0, 0, 600*sizemod, 600*sizemod])
%% Tree subsetting
nMu = 20;

[P, fatigue, arborLengths, m_csa]  = define_motor_pool(nMu, forceRange, fatigueRange, contractionRange);

n_arborLengths = ceil(arborLengths * 0.15 * n_nodes);
clf
mu_tree = select_MU_nodes(n_arborLengths, neuropoints);
% mu_tree = simple_MU_nodes(nMu, neuropoints);

for i = 1:nMu
    create_nRv(neuropoints, mu_tree{i});
    plot_tree(neuropoints, mu_tree{i} ,scale_gradient(i,:));
end
% plot_tree(neuropoints,

[x0, y0, r, h] = find_cyl(neuropoints);
sl = h/2;%r*1.2;
axis(axis_vals)
view(-45, 25)
%%
% export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\SubTree.png' -transparent -m3
spin_plot(200,[-sl sl -sl sl 0 h],'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\spin_tree')
% viscircles([x0,y0],r)

%% Treecutting
npnew = subset_tree(neuropoints, unique([mu_tree{:}]));
% for mu = 4:nMu
%     clf
%     plot_tree(npnew, scale_discrete(9,:), []);
%     plot_tree(npnew, scale_gradient(mod(mu,length(scale_gradient))+1,:), mu_tree{mu});
%     axis([-sl sl -sl sl 0 h])
%     view(-55, 15)
%     axis vis3d
%     pause;
% end

cycle_mu_plot(npnew,mu_tree,axis_vals,'\mu_cycle')
% spin_plot(100, [-sl sl -sl sl 0 h], [])
% export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\Ding2002Force.png' -transparent -m3
%% Place electrode in center and check activation of motor units
% mu = 6;
% tic
% r_pos = 0:1:r;
% steps = length(r_pos);
% np_activated = {};
% for dx = 1:steps
%     active_bool = zeros(nMu,1);
%     p_elec = [x0+r_pos(dx), y0, h/2];
%     
%     for mu = 1:nMu
%         active_bool(dx, mu) = is_activated(mu_tree{mu}, npnew, p_elec);
%     end
%     %     [mu_pts,dtrs] = dot_tree(npnew);
%     %     vxt = calc_Vext(p_elec/100, mu_pts/100, 10);
%     %     adj_mat = create_adj_matrix(dtrs, npnew);
%     %     D = 10e-6;
%     %     t_pulse = 0.5e-3;
%     %     activity = neuron_response(vxt'*1e-3,adj_mat,D,t_pulse);
%     
%     np_activated{dx} = unique([mu_tree{find(active_bool(dx,:))}]);
%     
%     figure(10)
%     clf
%     hold on
%     grid on
% %     
%     plot_tree(npnew, scale_discrete(4,:), []);
%     plot_tree(npnew, scale_discrete(2,:), [np_activated{dx}]);
% %     
% %     %     pause(3)
% %     %     plot3(mu_pts(find(activity),1),mu_pts(find(activity),2),mu_pts(find(activity),3),'.','Color', reed)
% %     %     plot3(mu_pts(~activity,1),mu_pts(~activity,2),mu_pts(~activity,3),'.','Color', bloo)
% %     plot3(p_elec(1),p_elec(2),p_elec(3),'o','Color', scale_discrete(8,:))
% %     plot3(x0+r_pos,repmat(p_elec(2),1,steps),repmat(p_elec(3),1,steps),'o','Color', scale_discrete(9,:))
% %     axis([-sl sl -sl sl 0 h])
% %     view(10, 11)
%     percent_force_recruitment(dx) = sum(P(find(active_bool(dx,:))))/sum(P)
% %     pause
% end
% toc

% r_pos = 4.5;
% 
% r_pos = 0:1:r;
% steps = length(r_pos);
% np_activated_thta = {};
% 
% for dx = 1:steps
%     active_bool = zeros(nMu,1);
%     p_elec = [x0+r_pos(dx), y0, h/2];
%     
%     for mu = 1:nMu
%         active_bool(dx, mu) = is_activated(mu_tree{mu}, npnew, p_elec);
%     end
%     %     [mu_pts,dtrs] = dot_tree(npnew);
%     %     vxt = calc_Vext(p_elec/100, mu_pts/100, 10);
%     %     adj_mat = create_adj_matrix(dtrs, npnew);
%     %     D = 10e-6;
%     %     t_pulse = 0.5e-3;
%     %     activity = neuron_response(vxt'*1e-3,adj_mat,D,t_pulse);
%     
%     np_activated{dx} = unique([mu_tree{find(active_bool(dx,:))}]);
%     
% %     figure(10)
% %     clf
% %     hold on
% %     grid on
% %     
% %     plot_tree(npnew, scale_discrete(1,:), []);
% %     plot_tree(npnew, scale_discrete(2,:), np_activated);
% %     
% %     %     pause(3)
% %     %     plot3(mu_pts(find(activity),1),mu_pts(find(activity),2),mu_pts(find(activity),3),'.','Color', reed)
% %     %     plot3(mu_pts(~activity,1),mu_pts(~activity,2),mu_pts(~activity,3),'.','Color', bloo)
% %     plot3(p_elec(1),p_elec(2),p_elec(3),'o','Color', scale_discrete(8,:))
% %     plot3(x0+r_pos,repmat(p_elec(2),1,steps),repmat(p_elec(3),1,steps),'o','Color', scale_discrete(9,:))
% %     axis([-sl sl -sl sl 0 h])
% %     view(10, 11)
%     percent_force_recruitment(dx) = sum(P(find(active_bool(dx,:))))/sum(P)
% %     pause
% end
% toc
% surf(activity)
% colormap summer
% shading interp
% light
%% NRV sequencing
% for i = 1:length(adj_mat)
%     clf
%     view(-45, 25)
%     axis([-sl sl -sl sl 0 h])
%     pnds = find(adj_mat(i,:));
%     plot_tree(npnew, scale_discrete(1,:),mu_tree{mu});
%     scatter3(mu_pts(pnds,1),mu_pts(pnds,2),mu_pts(pnds,3))
%     scatter3(mu_pts(i,1),mu_pts(i,2),mu_pts(i,3),'r')
%     pause
% end
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
%% Electrode Distances
% ~15 minutes
%Points in motor unit 1
figure(9)
clf
hold on
% tree1pts = dot_tree(neuropoints);
% tree1pts = dot_tree(npnew);
% plot_tree(tree1pts, bloo, [])
% dot_tree(neuropoints,mu_tree{2})

%meshgrid setup
THETA = linspace(0, 2*pi, 9);
THETA = THETA(1:end-1);
RHO = linspace(1, r*.7, 5);
Z = linspace(0,h,5);

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

%calculate distances from each cartesian point to each point in th tree

% dsts = pdist2(e_pts,tree1pts(:,1:3));
percent_force_recruitment = zeros(length(elec_x),1);
active_mu_bool = logical(zeros(length(elec_x),nMu));
np_mesh_activated = cell(nMu,1);
tic
for i = 1:length(elec_x)
    
    p_elec = [elec_x(i), elec_y(i), elec_z(i)];
    
    for mu = 1:nMu
        active_mu_bool(i, mu) = is_activated(mu_tree{mu}, npnew, p_elec);
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
    %     plot_tree(npnew, scale_discrete(1,:), []);
    %     plot_tree(npnew, scale_discrete(2,:), np_activated);
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
% plot_tree(neuropoints, bloo, [])
%get groups by rho and z
groupd = findgroups(rho,z);
%apply mean function to distances by groups of rho-z combinations
gmeans = splitapply(@mean,percent_force_recruitment,groupd);
distnc = reshape(gmeans, length(Z), length(RHO));

[rr, elec_z] = meshgrid(RHO,Z);
scatter3(rr(:),zeros(length(rr(:)),1), elec_z(:), distnc(:)*40+1,'r')
hold offc
view([0 0])
%%
toc
cycle_activation_plot2(npnew,mu_tree,active_mu_bool,e_pts, axis_vals,percent_force_recruitment','CycleActivations2.0')
%%
distnc(1) = 1;
figure(10)
clf
h = heatmap(RHO,Z, flip(distnc,1),'CellLabelColor','k', 'CellLabelFormat','%0.2g');
h.Title = 'Mean Proportion of Total Force Activated';
h.XLabel = 'Normalized Radial Position';
h.YLabel = 'Normalized Proximal-Distal Position';

