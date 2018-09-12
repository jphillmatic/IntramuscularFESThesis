%% Color Scales
global bloo
global reed
bloo = [5,112,176]./255;
reed = [215,48,31]./255;
scale_gradient0=[158,1,66
    213,62,79
    244,109,67
    253,174,97
    254,224,139
    171,221,164
    102,194,165
    50,136,189
    94,79,162]/255;

global scale_continuous
scale_continuous = flip([215,25,28
253,174,97
171,221,164
43,131,186]/256,1);

global scale_gradient
scale_gradient = [];
for i = 1:length(scale_gradient0)-1
    ncolors = linspaceNDim(scale_gradient0(i,:), scale_gradient0(i+1,:), max(3,ceil(nMu / (length(scale_gradient0)-1))))';
    scale_gradient = vertcat(scale_gradient,ncolors);
end
colorinds = round(linspace(1,length(scale_gradient),nMu));
scale_gradient = flip(scale_gradient(colorinds,:),1);
clear scale_gradient0

global scale_discrete
scale_discrete = [55,126,184
228,26,28
255,127,0
77,175,74
152,78,163
255,255,51
166,86,40
247,129,191
153,153,153]./255;