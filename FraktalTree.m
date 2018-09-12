function [neuropoints] = FraktalTree(n,prevNode,level,node) 
%% Input Arguments:
% n = # of maximum layers to iterate through
% a0 = starting position
% level = current level
% node = the node
if nargin == 1
    prevNode = 1;
    level = 1;
    node = 1;
end

if level == 1
    a0 = [0,0,0];
    text(a0(1), a0(2), a0(3), num2str(node))
    n0 = [0,0,0,0,0,1,0];
    parent = 1;
    clf;
    olddir = 0;
else
    n0 = [];
    a0 = prevNode(1:3);
    parent = prevNode(6);
    olddir = atan2(prevNode(2),prevNode(1));
end


% branch_area = 0.3 ^ (level - 1);
branch_area =  0;
rbranch = rand(1) + branch_area;

%87 axons exhibited 979 instances of branching ~~ 11.25 branches per axon
%979 instances of axonal branching, most (88.5%) were binary, 
%with progressively smaller fractions of higher degree branching 
%(tri-furcations 10.7%, 4-furcations 0.6%, 5-furcations 0.2%, see Figure S2)
straight_odds = 0.7;
nmj_odds = 0.2*(level/n)^0.5;
if level < 3
    nmj_odds = 0;
end
if rbranch < nmj_odds
    n_branches = 0;
elseif rbranch <= straight_odds + nmj_odds
    n_branches = 1;
else
    rnbranch = rand(1);
    if rnbranch <= 0.002 %0.2% of branches are quint
        n_branches = 5;
    elseif rnbranch <= 0.006 + 0.002 %.6% of branches are quad
        n_branches = 4;
    elseif rnbranch <= 0.107 + 0.006 + 0.002 %10.7% of branches are tri
        n_branches = 3;
    else %88.5% of branches are bi
        n_branches = 2; 
    end
end

if level == 1
    n_branches = 3 + ceil(rand(1)*3);
end

%Branch Length
l_branch = randn(n_branches,1) + 4;

%Phi = radial angle
%Chi = vertical angle
chi = pi/2 - abs(pi/6*rand(1));
phi = olddir + pi/2*randn(1);
if n_branches > 1
    phi0 = mod(phi, 2*pi/n_branches);
    for i = 2:n_branches
        chi(i) = pi/4 + pi / 8 * randn(1);
        phi(i) = olddir + phi0 + pi/4*randn(1)/n_branches;
    end
end

if level == 1
    chi = pi / 5 + pi/12*randn(n_branches,1);
    chi(end) = pi / 2 + pi/12*randn(1);
    phi = [linspace(0, 2*pi, n_branches) + pi/4/n_branches*randn(1,n_branches)]';
    l_branch(end) = 4;
end

% chi = pi/2 - abs(pi/4/level^0.5  + pi / 8 * randn(n_branches,1));

% if uniform_tree
%     n_branches = 3;
%     l_branch = ones(n_branches,1);
%     chi = pi/2.5 * ones(n_branches,1);
%     phi = linspace(0, 2 * pi, n_branches+1);
%     phi = phi(1:end-1);
% end

newpoints = [];
nextpoints = [];
nxp = [];

for i = 1:n_branches
    %Generate New x,y,z values
    %Phi = radial angle
    %Chi = vertical angle
    xa = l_branch(i) * cos(phi(i)) * cos(chi(i));
    ya = l_branch(i) * sin(phi(i)) * cos(chi(i));
    za = l_branch(i) * sin(chi(i));
    a1 = [xa,ya,za] + a0;
    %Increment Node#
    if ((i == 1) || (isempty(nextpoints)))
        node = node + 1;
    else
        node = max(max(nextpoints(:,6)) , node)+1;
    end
    
    %Contain new node data in np
    np = horzcat(a1, level, i, node, parent);
    
    %Add np to the list of points at this level
    newpoints = vertcat(newpoints, np);
    
    %Draw lines from original point to each new point
%     line([a0(1), a1(1)], [a0(2), a1(2)], [a0(3), a1(3)],'Color', [.75 .25 .25])
%     text(a1(1), a1(2), a1(3), num2str(node))
    %If not at the desired level, call function on all nodes
    if level < n
        %Call FraktalTee using current node as origin
        %Increment level and node by 1 in function call
        nxp  = FraktalTree(n,np,level+1,node);
        
        %Create list of nodes from furure levels
        nextpoints = vertcat(nextpoints,nxp);
    end  
end

%When maximum level has been reached, return newpoints, and future points
neuropoints = vertcat(n0, newpoints, nextpoints);