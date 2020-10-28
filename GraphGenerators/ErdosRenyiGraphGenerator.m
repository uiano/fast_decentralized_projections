classdef ErdosRenyiGraphGenerator < GraphGenerator
	
	properties % required by parent classes
		c_parsToPrint  = {'ch_name','prob_edge','num_nodes'};
		c_stringToPrint  = {'','Edge prob.',''};
		c_patternToPrint = {'%s%s random graph','%s = %g','%s%d vertices'};
	end
	
	properties(Constant)
		ch_name = 'Erdos-Renyi';
	end
	
	properties
		prob_edge;
		num_nodes;
		prob_selfLoop = 0;
		b_symmetric;
	end
	
	methods
		
		function obj = ErdosRenyiGraphGenerator(varargin)
			% Constructor
			obj@GraphGenerator(varargin{:});
			
		end
		
		function graph = realization(obj)
			% Output:
			% GRAPH       Object of class Graph which contains a
			%             realization of a symmetric Erdos-Renyi random graph
			%
			
			assert(~isempty(obj.prob_edge));
			assert(~isempty(obj.num_nodes));
			
			m_adjacency = rand(obj.num_nodes) < obj.prob_edge;
			if obj.b_symmetric==0
				for ind_diag=1:obj.num_nodes
					m_adjacency(ind_diag,ind_diag)=1;
				end
			end
			if obj.b_symmetric==1
				m_adjacency = m_adjacency - diag(diag(m_adjacency));
				m_adjacency = triu(m_adjacency) + triu(m_adjacency)';
				
				if obj.prob_selfLoop
					v_diagonal = rand(1, obj.num_nodes)<=obj.prob_selfLoop;
					m_adjacency = m_adjacency + diag(v_diagonal);
				end
			end
			graph = Graph('m_adjacency',m_adjacency);
			if obj.b_checkConnected==1
				[	graph] = graph.randomlyConnect();
			end
		end
	end
end

