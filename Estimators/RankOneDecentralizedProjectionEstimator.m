classdef RankOneDecentralizedProjectionEstimator < DecentralizedProjectionEstimator
	
	properties
		
		b_diag %determines the diagonal elements of the shift matrix 
		b_nodeVariant = 1; % set to 0 to use node-invariant filters
	end
	
	
	methods		
		
		function [m_Shift,m_projection]=findShiftMatrix(obj,m_basisSubspace,graph)
					
			% INPUT:
			%   m_basisSubspace: full-rank num_nodes x 1 matrix
			%       whose columns span the signal subspace
			%   graph: object of class Graph
			% OUTPUT:
			%   m_shift : num_nodes x num_nodes 
			%       matrix where m_shift
			%        is the shift matrix
			
			% input processing
			assert(size(m_basisSubspace,2)==1);%The dimension of the subspace must be 1
			num_nodes=size(m_basisSubspace,1);
			m_projection=m_basisSubspace*m_basisSubspace';
			S=zeros(num_nodes,num_nodes);
			% obtaining the connected edges
			m_connected=graph.getConnectedEdges();
			% find the shift matrix
			for ind=1:size(m_connected,2)
				if m_connected(1,ind)~=m_connected(2,ind)
					S(m_connected(1,ind),m_connected(2,ind))=-m_basisSubspace(m_connected(1,ind))*m_basisSubspace(m_connected(2,ind));
				else
					[s_row,s_col]=find(m_connected(1,:)==m_connected(1,ind));
					s_sum=0;
					for ind_col=1:size(s_col,2)
						s_loc=m_connected(2,s_col(ind_col));
						if m_connected(2,s_loc)~=m_connected(1,ind)
							s_sum=s_sum+m_basisSubspace(s_loc)*m_basisSubspace(s_loc);
						end
						S(m_connected(1,ind),m_connected(2,ind))=obj.b_diag+s_sum;
					end
				end
			end
			m_Shift=S;%the shift matrix
		end
	end
end