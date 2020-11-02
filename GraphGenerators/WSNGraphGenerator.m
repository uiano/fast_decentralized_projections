classdef  WSNGraphGenerator < GraphGenerator
	
	properties
		num_nodes;% the number of nodes
		threshold;% If the distance between each nodes is lower than this parameter, we conclude that an edge exists between nodes
	end
	
	methods
		
		function obj = WSNGraphGenerator(varargin)
			% Constructor
			obj@GraphGenerator(varargin{:});
			
		end
		
		function graph = realization(obj)
			
			% it is a graph generator for WSNs. Firstly, it produces a number
			% of nodes locations(equals num_nodes) randomly in an area
			% (uniform over a 1 by 1 square). Then it compares the distance
			% between each nodes. If that distance is lower than a
			% threshold, it concludes that those nodes have a connection with each other.
			% Finally, by this algorithm, we can find thhe adjacency matrix.
			%
			% OUTPUT:
			%   graph :the adjacency matrix
			
			%input processing
			x_coordinate=rand(obj.num_nodes,1);%x coordinate of nodes
			y_coordinate=rand(obj.num_nodes,1);%y coordinate of nodes
			m_adjacency=ones(obj.num_nodes,obj.num_nodes);
			
			%program
			for ind_row=1:obj.num_nodes
				for ind_column=ind_row:obj.num_nodes
					distance = sqrt((x_coordinate(ind_row) - x_coordinate(ind_column))^2 + (y_coordinate(ind_row) - y_coordinate(ind_column))^2);%distance between each pairs of nodes
					if distance >obj.threshold   %if it is correct then there is an edge between each pairs
						m_adjacency(ind_row,ind_column)=0;
						m_adjacency(ind_column,ind_row)=0;
					end
				end
			end
			graph = Graph('m_adjacency',m_adjacency);
			if obj.b_checkConnected==1
				[	graph] = graph.randomlyConnect();
			end
		end
	end
end