classdef NerveTreeClass
    properties
        Nodes
        Parents
        Positions
    end
    methods
        function [chld] = find_children(obj,node)
            chld = obj.Nodes(obj.Parents == node);
        end
    end
end
