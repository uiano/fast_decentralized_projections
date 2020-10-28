classdef GraphGenerator %< ProcessingBlock
	% Subclasses generate graphs either from real data or randomly
	
	properties(Constant)
	end
	
	properties
		b_checkConnected=1;%if b_connected=1, then a connected graph is generated
	end
	
	methods
		
		function obj = GraphGenerator(varargin)
			%obj@ProcessingBlock(varargin{:});
		end
		
	end
	
	methods(Abstract)
		
		graph = realization(obj);
		% GRAPH      Object of class Graph
		
	end
	
end

