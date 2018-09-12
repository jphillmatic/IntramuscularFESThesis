[mu_thru] = findPathWeights(npnew, mu_tree);
mu_thru(:,2) = (mu_thru(:,2) / max(mu_thru(:,2))).^0.4;
figure(1)
clf
subplot(2,3,[1,2,4,5])
for i = 1:nMu
    pmu = nMu + 1 - i ;
%     create_nRv(neuropoints, mu_tree{i});
    plot_tree(npnew, mu_tree{pmu},scale_gradient(i,:), mu_thru);
end
title('All Motor Units')
view(-15, 15)
axis tight
axis([-sl sl -sl sl 0 ht])

[mu_thru] = findPathWeights(npnew, mu_tree);
mu_thru(:,2) = (mu_thru(:,2) / max(mu_thru(:,2))).^0.5;

subplot(2,3,[3])
plot_tree(npnew, [], scale_discrete(9,:), mu_thru);
plot_tree(npnew, mu_tree{1}, scale_gradient(1,:), mu_thru);
title('Smallest Motor Unit')
axis tight
axis([-sl sl -sl sl 0 ht])

view(-105, 15)

subplot(2,3,[6])
plot_tree(npnew, [], scale_discrete(9,:), mu_thru);
plot_tree(npnew, mu_tree{nMu}, scale_gradient(nMu,:), mu_thru);
title('Largest Motor Unit')
axis tight
axis([-sl sl -sl sl 0 ht])
view(75, 15)


whitebg('white')
set(gcf, 'Color', [1 1 1])
% whitebg('white')
% set(gcf, 'Color', [0 0 0])
set(gcf, 'Position', [0, 0, 900*1.5, 600*1.5])
[x0, y0, r, ht] = find_cyl(npnew);
sl = ht/4;
axis_vals = [-sl sl -sl sl 0 ht];
axis(axis_vals)

% export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\individualmus.png' -transparent -m3
