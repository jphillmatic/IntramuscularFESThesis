soleus_length = 39.3; %cm
soleus_volume = 356+58; %cm3
soleus_radius = sqrt(soleus_volume / soleus_length / pi); %cm

for runx = 1:25
    nMu = 40;
    nNMJ = nMu *4;
    
    n = 5;
    nminorbranches = 0;
    
    unbalanced = true;
    while (nminorbranches < 18) || unbalanced1 || unbalanced2
        np0 = FraktalTwigTree872018(20);
        [x0, y0, r, ht] = find_cyl(np0);
        
        
        np0 = sortrows(np0,6);
        np0(:,3) = np0(:,3);
        ht_mod = soleus_length/ht;
        np0(:,3) = np0(:,3)*ht_mod;
        
        [tta rda] = cart2pol(np0(:,1)+x0, np0(:,2)+y0);
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
        znmj = rand(nNMJ,1)*3*ht/10 + ht/5;
        
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
    mu_thru(:,2) = (mu_thru(:,2) / max(mu_thru(:,2))).^0.6;

    for mu = 1:nMu
        npmu = subset_tree(neuropoints,mu_tree{mu});
        for i = 1:length(npmu(:,6))
            nd(i,mu) = length(find_daughters(npmu(i,6), npmu));
        end
    end
    
    tot_brnchs = sum(sum(nd>1));
    brnchs_rate = sum(sum(nd>1))/nMu;
    tot_xb(runs,:) = [sum(sum(nd==2))/tot_brnchs, sum(sum(nd==3))/tot_brnchs, sum(sum(nd==4))/tot_brnchs, sum(sum(nd==5))/tot_brnchs];
    
    % ~15 minutes
    %Points in motor unit 1

    % tree1pts = dot_tree(neuropoints);
    % tree1pts = dot_tree(npnew);
    % plot_tree(tree1pts, [], bloo)
    % dot_tree(neuropoints,mu_tree{2})
    
    %meshgrid setup
    THETA = linspace(0, 2*pi, 9);
    THETA = THETA(1:end-1);
    RHO = linspace(r*0.1, r, 5);
    Z = linspace(0,ht*0.9,5);
   
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
            active_mu_bool(i, mu) = is_activated(mu_tree{mu}, npnew, p_elec);
        end
        percent_force_recruitment(runs,i) = sum(P(active_mu_bool(i,:)))/sum(P);
        
        np_mesh_activated{i} = unique([mu_tree{find(active_mu_bool(i,:))}]);

    end
    % [pt2treemin,ys] = min(dsts,[],2);
    
    % scatter3(e_pts(:,1), e_pts(:,2), e_pts(:,3), '.', 'MarkerEdgeAlpha',0.75)
    % plot_tree(neuropoints, [], bloo)
    %get groups by rho and z
    groupd = findgroups(rho,z);
    %apply mean function to distances by groups of rho-z combinations
    gmeans = splitapply(@mean,percent_force_recruitment(runs,:),groupd);
    distnc = flip(reshape(gmeans, length(Z), length(RHO)),1);
    
    distncruns(:,:,runs) = distnc;
    [rr, elec_z] = meshgrid(RHO,Z);

end
