function spin_plot(n_frames,axis_vals,fname)

% set(gca,'nextplot','replacechildren','visible','off')
% set(gcf, 'Renderer', 'zbuffer')
% set(gcf, 'Color', [1 1 1])

xt0 = xticks(gca);
zt0 = zticks(gca);

f = getframe(gcf);
% direction = [0 0 1];
[im,map] = rgb2ind(f.cdata,256,'nodither');
axis vis3d
im(1,1,1,20) = 0;
i = 1;
az =20;
% ax = gcf;
for thta = linspace(-45, 360-45, n_frames)
    %   az = 5*cos(thta/2/pi/6)+15;
    view([thta, az])
    
    axis(axis_vals)

    xticks(xt0)
    yticks(xt0)
    zticks(zt0)
    
    f = getframe(gcf);
    im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    i = i + 1;
end

if ~isempty(fname)
    imwrite(im,map,[fname '.gif'],'DelayTime',0,'LoopCount',inf)
end