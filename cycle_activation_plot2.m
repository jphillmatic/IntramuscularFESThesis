function cycle_activation_plot2(neuropoints,mu_tree,mu_bool,e_pts,axis_vals,cscale,fname)

global scale_discrete
global scale_continuous
global scale_gradient
global bloo
nMu = size(mu_tree,2);
nEpos = size(e_pts,1);
figure(1)
clf
si_top = [0 0 ; -90, 90];
% for i = 1:2
%     subplot(1,2,i)
%     plot_tree(npnew, scale_gradient(1,:), []);
%     hold on
%     scatter3(e_pts(:,1), e_pts(:,2), e_pts(:,3), '.', 'MarkerEdgeAlpha',0.75)
%     scatter3(e_pts(i,1), e_pts(i,2), e_pts(i,3), 1000, scale_discrete(4,:), '.' )
%
%     scatter3(e_pts(1,1), e_pts(1,2), e_pts(1,3), 50, scale_discrete(2,:))
%     for j = 1:length(scale_continuous)
%         plot_tree(npnew,scale_continuous(i,:), mu_tree{1});
%     end
%
%     view(si_top(i,:))
%
%     whitebg('black')
% end
% axis([axis_vals])
% f = getframe(gcf);
%
% set(gcf, 'Color', [0 0 0])
% set(gcf, 'Position', [0, 0, 1000*1.5, 400*1.5])

[mu_thru] = findPathWeights(neuropoints, mu_tree);
mu_thru(:,2) = (mu_thru(:,2) / max(mu_thru(:,2))).^0.5;

% im(1,1,1,nEpos) = 0;

for i = 1:nEpos
    clf
    active_nodes = unique([mu_tree{mu_bool(i,:)}]);
    
    for j = 1:2
        subplot(1,2,j)
        
       plot_tree(neuropoints, [], scale_discrete(9,:), mu_thru);
        hold on
        
        if ~isempty(active_nodes)
            %         plot_tree(npnew, scale_continuous(ceil(cscale(i)*length(scale_continuous)),:), active_nodes);
            plot_tree(neuropoints, active_nodes, scale_discrete(2,:), mu_thru);
        end
        scatter3(e_pts(:,1), e_pts(:,2), e_pts(:,3),10, bloo, 'filled');
        scatter3(e_pts(i,1), e_pts(i,2), e_pts(i,3),75, scale_discrete(4,:), 'filled');
        
        
%         plot3(e_pts(:,1), e_pts(:,2), e_pts(:,3), '.')
%         plot3(e_pts(i,1), e_pts(i,2), e_pts(i,3), '.', 1000, 'color', scale_discrete(4,:) )
        
        
        
%         legend("Nerve Tree", "Activated Branches", "Electrode Position",'Location','NorthEastOutside')        
        view(si_top(j,:))
%         axis([axis_vals])
    end
    
    whitebg('black')
    set(gcf, 'Color', [0 0 0])
    set(gcf, 'Position', [0, 0, 1000*1.5, 400*1.5])
    
    f = getframe(gcf);

    
    pause(.5)
    
    if i == 1
        [im,map] = rgb2ind(f.cdata,256,'nodither');
    else
        im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    end
    
    
    %     clf
    %     subplot(1,2,1)
    %     plot_tree(npnew, scale_gradient(1,:), []);
    %     hold on
    %     active_nodes = unique([mu_tree{find(mu_bool(i,:))}]);
    %     if ~isempty(active_nodes)
    % %         plot_tree(npnew, scale_continuous(ceil(cscale(i)*length(scale_continuous)),:), active_nodes);
    %         plot_tree(npnew, scale_discrete(2,:), active_nodes);
    %
    %     end
    %
    %     scatter3(e_pts(:,1), e_pts(:,2), e_pts(:,3), '.', 'MarkerEdgeAlpha',0.75)
    %     scatter3(e_pts(i,1), e_pts(i,2), e_pts(i,3), 1000, scale_discrete(4,:), '.' )
    %     axis([axis_vals])
    %     set(gcf, 'Color', 'black')
    %     whitebg('black')
    %     legend("Nerve Tree", "Activated Nodes", "Electrode Position")
    %
    %     view([-90, 88])
    %
    %     subplot(2,2,1)
    %     plot_tree(npnew, scale_gradient(1,:), []);
    %     hold on
    %     active_nodes = unique([mu_tree{find(mu_bool(i,:))}]);
    %     if ~isempty(active_nodes)
    % %         plot_tree(npnew, scale_continuous(ceil(cscale(i)*length(scale_continuous)),:), active_nodes);
    %         plot_tree(npnew, scale_discrete(2,:), active_nodes);
    %
    %     end
    %
    %     scatter3(e_pts(:,1), e_pts(:,2), e_pts(:,3), '.', 'MarkerEdgeAlpha',0.75)
    %     scatter3(e_pts(i,1), e_pts(i,2), e_pts(i,3), 1000, scale_discrete(4,:), '.' )
    %     axis([axis_vals])
    %
    %     set(gcf, 'Color', [0 0 0])
    %     whitebg('black')
    %     legend("Nerve Tree", "Activated Nodes", "Electrode Position")
    %     view([0, 0])
    
    %     set(gcf, 'Position', [0, 0, 1000*1.5, 400*1.5])
    
    %     f = getframe(gcf);
    %     im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    %     pause;
    
end

if ~isempty(fname)
    imwrite(im,map,['C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Project Files\Jake Model\Figures\' fname '.gif'],'DelayTime',0.25,'LoopCount',inf)
end