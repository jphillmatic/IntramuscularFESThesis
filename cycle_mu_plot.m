function cycle_mu_plot(npnew,mu_tree,axis_vals,fname)

% set(gca,'nextplot','replacechildren','visible','off')
% set(gcf, 'Renderer', 'zbuffer')
% set(gcf, 'Color', [1 1 1])
global scale_discrete
global scale_gradient
nMu = size(mu_tree,2);
figure(1)
clf
az =20;

for mu = 1:nMu
    plot_tree(npnew,scale_gradient(mod(mu,length(scale_discrete))+1,:), mu_tree{mu});
end
axis([axis_vals])
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');

im(1,1,1,nMu) = 0;

    i = 1;
ax = gcf;
for mu = 1:nMu

    plot_tree(npnew, scale_discrete(9,:), []);
    hold on
    plot_tree(npnew, scale_gradient(mod(mu,length(scale_discrete))+1,:), mu_tree{mu});
    axis([axis_vals])
    view([0, 90])
    
    f = getframe(gcf);
    im(:,:,1,mu) = rgb2ind(f.cdata,map,'nodither');
    clf
end

if ~isempty(fname)
    imwrite(im,map,['C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Project Files\Jake Model\Figures\' fname '.gif'],'DelayTime',0.3,'LoopCount',inf)
end