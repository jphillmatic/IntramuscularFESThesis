function [mu_bool] = cycle_current_plot(neuropoints,mu_tree,i_stims,fname)


global scale_discrete

nMu = size(mu_tree,2);
nIstims = length(i_stims);
figure(1)
clf
ysps = ceil(nIstims/4);
subplot(ysps,ceil(nIstims/ysps),1)

[x0, y0, r, ht] = find_cyl(neuropoints);
[p_elec(1),p_elec(2), p_elec(3)]  = pol2cart(pi/4, r*2/3,ht/2);
% p_elec = [2*r/3, 0, ht/2]

for i = 1:nIstims    
    for mu = 1:nMu
        mu_bool(i, mu) = is_activated(mu_tree{mu}, neuropoints, p_elec, 0.5e-3, i_stims(i));
    end
end

[mu_thru] = findPathWeights(neuropoints, mu_tree);
mu_thru(1:2,2) = mu_thru(1:2,2)/2;
mu_thru(:,2) = (mu_thru(:,2) / max(mu_thru(:,2))).^0.5;

scatter3(2*r/3, 0, ht/2, 10, scale_discrete(4,:), 'filled' )
plot_tree(neuropoints, [], scale_discrete(9,:)*1.1);
plot_tree(neuropoints, [], scale_discrete(2,:), mu_thru)

% legend("Electrode Position", "Nerve Tree", "Activated Nodes")

whitebg('white')
set(gcf, 'Renderer', 'opengl')
set(gcf, 'Color', [1 1 1])
% set(gcf, 'Position', [0, 0, 600*1.5, 600*1.5])
set(gcf, 'Position', [0, 0, 1600, 900])



sl = ht/4;
axis_vals = [-sl sl -sl sl 0 ht];
axis([axis_vals])
view(0, 0)


f = getframe(gcf);

% [im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,nIstims) = 0;

i = 1;
ax = gcf;
for i = 1:nIstims
%     clf
    subplot(ysps, ceil(nIstims/ysps),i)
    
    scatter3(p_elec(1), p_elec(2), p_elec(3), 5000*i_stims(i), scale_discrete(4,:), 'filled' )
    
    active_nodes = unique([mu_tree{find(mu_bool(i,:))}]);
    
    plot_tree(neuropoints, [], scale_discrete(9,:)*1.1, mu_thru);
    hold on
    
    if ~isempty(active_nodes)
        %                 plot_tree(neuropoints, active_nodes, scale_continuous(ceil(cscale(i)*length(scale_continuous)),:), mu_thru);
        plot_tree(neuropoints, active_nodes, scale_discrete(2,:), mu_thru);
        
    end
    
    
    
%     legend("Electrode Position", "Nerve Tree", "Activated Nodes")
    
    title(['Motor Unit Activation, Stimulus Current = ' [num2str(i_stims(i)*1e3)] 'mA'])
    
    axis([axis_vals])
    view(0, 0)
    
    
    
%     f = getframe(gcf);
%     im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    %     pause;
    
end

if ~isempty(fname)
    export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\stim_current_plot.png' -transparent -m3
end
% if ~isempty(fname)
%     imwrite(im,map,['C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Project Files\Jake Model\Figures\' fname '.gif'],'DelayTime',1,'LoopCount',inf)
% end