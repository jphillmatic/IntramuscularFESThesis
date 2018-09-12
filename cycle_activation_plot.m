function cycle_position_plot(neuropoints,mu_tree,mu_bool,e_pts,axis_vals,cscale,fname)

global scale_discrete

nMu = size(mu_tree,2);
nEpos = length(e_pts);
figure(1)
clf

[mu_thru] = findPathWeights(neuropoints, mu_tree);
mu_thru(:,2) = (mu_thru(:,2) / max(mu_thru(:,2))).^0.5;

plot_tree(neuropoints, [], scale_discrete(9,:)*1.1);
plot_tree(neuropoints, [], scale_discrete(2,:), mu_thru)
scatter3(e_pts(:,1), e_pts(:,2), e_pts(:,3), 10, scale_discrete(1,:), 'filled')
scatter3(e_pts(1,1), e_pts(1,2), e_pts(1,3), 100, scale_discrete(4,:), 'filled' )

legend("Nerve Tree", "Activated Nodes", "Electrode Position", "Possible Positions")

whitebg('black')
set(gcf, 'Renderer', 'opengl')
set(gcf, 'Color', [0 0 0])
set(gcf, 'Position', [0, 0, 600*1.5, 600*1.5])

[x0, y0, r, ht] = find_cyl(neuropoints);
sl = ht/4;
axis_vals = [-sl sl -sl sl 0 ht];
axis([axis_vals])
view(10, 20)


f = getframe(gcf);

[im,map] = rgb2ind(f.cdata,256,'nodither');

im(1,1,1,nEpos) = 0;

i = 1;
ax = gcf;
for i = 1:nEpos
    clf
    active_nodes = unique([mu_tree{find(mu_bool(i,:))}]);
    
    plot_tree(neuropoints, [], scale_discrete(9,:)*1.1, mu_thru);
    hold on
    
    if ~isempty(active_nodes)
        %plot_tree(neuropoints, active_nodes, scale_continuous(ceil(cscale(i)*length(scale_continuous)),:), mu_thru);
        plot_tree(neuropoints, active_nodes, scale_discrete(2,:), mu_thru);
        
    end
    
    scatter3(e_pts(:,1), e_pts(:,2), e_pts(:,3), 10, scale_discrete(1,:), 'filled')
    scatter3(e_pts(i,1), e_pts(i,2), e_pts(i,3), 100, scale_discrete(4,:), 'filled')
    
    title('Motor Unit Activation by Electrode Position')
    
    axis([axis_vals])
    view(10, 20)
    
    f = getframe(gcf);
    im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    %     pause;
    
end

if ~isempty(fname)
    imwrite(im,map,['C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Project Files\Jake Model\Figures\' fname '.gif'],'DelayTime',1,'LoopCount',inf)
end