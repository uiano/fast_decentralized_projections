classdef GossipingDecentralizedProjectionEstimator < DecentralizedProjectionEstimator
	
	properties
		
		b_verbose = 0; % set to 1 to print info on the screen while executing a method.
		
	end
	
	
	methods
		
		function [m_signalEstimates]=estimate(obj,noisySignal,m_basisSubspace,graph,num_localExchanges)
			
			%
			% It estimates the object signal by using the shift matrix (m_shift) obtained by the Gossiping method.
			%
			% INPUT:
			%   noisySignal : num_nodes x 1 vector with the observations
			%   m_basisSubspace: full-rank num_nodes x dim_subspace matrix
			%       whose columns span the signal subspace
			%   graph : object of class Graph
			%   num_localExchanges : number of times a node shares its
			%       estimate with its neighbors.
			% OUTPUT:
			%   m_signalEstimates : num_nodes x num_local_exchanges matrix
			%       where the n-th column is the signal estimate
			%       after n local exchanges.
			
			% input processing
			num_nodes=size(m_basisSubspace,1);
			
			%find The shift matrix
			[m_Shift]=obj.findShift(m_basisSubspace,graph);
			
			%find the signal estimation
			m_signalEstimates = NaN(num_nodes,num_localExchanges);
			m_signalEstimates(:,1) = m_Shift*noisySignal;
			for ind_localExchanges = 1:num_localExchanges-1
				m_signalEstimates(:,ind_localExchanges+1) = m_Shift*m_signalEstimates(:,ind_localExchanges);
			end
		end
		
		function [t_filterMatrices,m_Shift]=getFilterMatrices(obj,m_basisSubspace,graph,num_localExchanges)
			
			%
			% it Obtains a filter by using the gossiping method (m_Shift^(ind_localExchanges)).
			%
			% INPUT:
			%   m_basisSubspace: full-rank num_nodes x dim_subspace matrix
			%       whose columns span the signal subspace
			%   graph : object of class Graph
			%   num_localExchanges : number of times a node shares its
			%       estimate with its neighbors.
			% OUTPUT:
			%   t_filterMatrices : num_nodes x num_nodes x
			%       num_local_exchanges tensor where
			%       t_filterMatrices(:,:,n) is the filter matrix
			%       after n local exchanges.
			
			% input processing
			num_nodes=size(m_basisSubspace,1);
			
			
			%find The Filter
			t_filterMatrices=NaN(num_nodes,num_nodes,num_localExchanges);
			[m_Shift]=obj.findShift(m_basisSubspace,graph);
			for ind_localExchanges = 1:num_localExchanges
				t_filterMatrices(:,:,ind_localExchanges)=m_Shift^(ind_localExchanges);
			end
		end
		
		function b_feasible = isFeasible(obj,m_basisSubspace,graph)
			%
			% Returns 1 if the optimization solved to find the shift matrix
			% is feasible; 0 otherwise.
			%
			[~,b_feasible] = obj.findShift(m_basisSubspace,graph);
			
		end
		
		function [m_Shift,b_feasible]=findShift(obj,m_basisSubspace,graph)
			
			%
			% It obtains the shift matrix of the Gossiping method and checks whether the method (CVX) codes is feasible or not.
			%  if b_feasible=0 means the method is not feasible and it gets
			%  identity matrix as m_shift. If b_feasible=1 means the method
			%  is feasible and it returns m_shift via the Gossiping method.
			
			% INPUT:
			%   m_basisSubspace: full-rank num_nodes x dim_subspace matrix
			%       whose columns span the signal subspace
			%   graph: object of class Graph
			% OUTPUT:
			%   m_shift : num_nodes x num_nodes x
			%       matrix where m_shift
			%        is the shift matrix
			
			% input processing
			[~,U_perp] = obj.findOrthonormalBases(m_basisSubspace);
			num_nodes=size(U_perp,1);
			dim_orthogonalsubspace=size(U_perp,2);
			
			% obtaining the missing edges
			m_missingEdges=graph.getMissingEdges();
			num_zero=size(m_missingEdges,2);
			
			% find the shift matrix
			if obj.b_verbose
				cvx_begin sdp
			else
				cvx_begin sdp quiet
			end
			variable L_tilde(dim_orthogonalsubspace,dim_orthogonalsubspace) symmetric
			variable gama
			variable mu_tilde
			minimize(gama)
			subject to
			L_tilde-eye(dim_orthogonalsubspace) >=0
			L_tilde-gama*eye(dim_orthogonalsubspace) <=0
			mu_tilde==trace(L_tilde)
			mu_tilde>=0
			S=U_perp*L_tilde*U_perp';
			for ind_zeronodes=1:num_zero
				S(m_missingEdges(1,ind_zeronodes),m_missingEdges(2,ind_zeronodes))==0;
			end
			cvx_end
			
		%	 find optimal step size and weight matrix
			if strcmp(cvx_status,'Infeasible')==1
				b_feasible = 0;
				m_Shift=zeros(num_nodes);
			else
				b_feasible = 1;
				L_1=L_tilde*mu_tilde;
				eigen_L_1=sort(eig(L_1),'ascend');
				epsil=2/(eigen_L_1(1)+eigen_L_1(dim_orthogonalsubspace));
				L=mu_tilde*U_perp*L_tilde*U_perp';
				m_Shift=eye(num_nodes)-epsil*L;
			end
		end
		
	end
end