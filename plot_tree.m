function plot_tree(neuropoints, varargin)
% global scale_discrete
global bloo

figure(1)
hold on
grid on
xlabel('X (cm)')
ylabel('Y (cm)')
zlabel('Z (cm)')
title(' ')
sizemod = 1.5;
% set(gcf, 'Position', [0, 0, 600*sizemod, 600*sizemod])
% set(gca,'color',[1 1 1])
% rnd0 = rand(1)/30;


switch nargin
    case 1
        pclr = bloo;
        nodes = [];
        lweights = [];
    case 2
        nodes = varargin{1};
        pclr = bloo;
        lweights = [];
    case 3
        nodes = varargin{1};
        pclr = varargin{2};
        lweights = [];
    case 4
        nodes = varargin{1};
        pclr = varargin{2};
        lweights = varargin{3};
end

if ~isempty(nodes)
    neuropoints = subset_tree(neuropoints, nodes);
else
    nodes = neuropoints(:,6);
end

thiccnes = 6;
dotmod = 6;

scatter3(0, 0, 0,dotmod*thiccnes, pclr, 'filled')

if isempty(lweights)
    scatter3(neuropoints(:,1), neuropoints(:,2), neuropoints(:,3),dotmod*thiccnes./neuropoints(:,4).^2, pclr, 'filled')
else
    scatter3(neuropoints(:,1), neuropoints(:,2), neuropoints(:,3),thiccnes*dotmod*lweights(ismember(lweights(:,1),nodes),2), pclr, 'filled')
end
% text(neuropoints(:,1), neuropoints(:,2), neuropoints(:,3), num2str(neuropoints(:,6)))

parents = unique(neuropoints(:,7));

for i = 1:length(parents)
    p = parents(i);
    
    pnp = neuropoints(neuropoints(:,6) == p,:);
    
    daughters = find_daughters(p, neuropoints);
    
    if ~(isempty(daughters) || isempty(pnp))
        for j = 1:length(daughters)
            d = daughters(j);
            dnp = neuropoints(neuropoints(:,6) == d,1:4);
            
            if nargin < 4
                line([pnp(1), dnp(1)], [pnp(2), dnp(2)], [pnp(3), dnp(3)],'LineWidth',thiccnes/dnp(4)^0.75, 'Color', pclr)
            else
                line([pnp(1), dnp(1)], [pnp(2), dnp(2)], [pnp(3), dnp(3)],'LineWidth',1.25*thiccnes * lweights(lweights(:,1)==d,2), 'Color', pclr)
            end
           
        end
    end
end
