classdef DecentralizedProjectionEstimator < matlab.mixin.Heterogeneous
	%DECENTRALIZEDPROJECTIONESTIMATOR Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		b_second=0 % determines that for the succeesive method, the second GF
		% is used or not
	end
	
	
	methods
		
		function b_feasible = isFeasible(obj)
			error('Method isFeasible() not implemented by subclass of DecentralizedProjectionEstimator');
			
		end
		
		function [ t_filterMatrices, m_shift] = getFilterMatrices(obj,m_basisSubspace,graph,num_localExchanges)
			%
			% Computation of graph filter matrices for decentralized subspace projection.
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
			%   t_normalizedError : 1 x 1 x num_local_exchanges tensor
			%       where the n-th entry is the error of t_filterMatrices(:,:,n)
			%       relative to the projection matrix
			%   m_shift : num_nodes x num_nodes shift matrix.
			num_nodes = size(m_basisSubspace,1);
			t_filterSecond=zeros(num_nodes,num_nodes,num_localExchanges);
			% Compute the shift matrix
			[m_shift,m_projection] = ...
				obj.findShiftMatrix(m_basisSubspace,graph);
			
			% Compute filter matrices
			t_filterMatrices = NaN(num_nodes,num_nodes,num_localExchanges);
			for ind_localExchanges = 1:num_localExchanges
				if norm(m_shift-eye(num_nodes))<1e-6
					t_filterMatrices(:,:,ind_localExchanges) = zeros(num_nodes); % good when r small.
				else
					[~,~,m_filter] = obj.findCoefficientsForOrder(m_shift,m_projection,min(ind_localExchanges,num_nodes-1));
					t_filterMatrices(:,:,ind_localExchanges) = m_filter;
				end
			end
			if obj.b_second==1
				[m_filterSecond]=obj.findSecondShift(m_projection-t_filterMatrices(:,:,end),graph);
				num_layer=size(m_filterSecond,3);
				for ind=1:num_layer
					t_filterSecond(:,:,end-num_layer+ind)=m_filterSecond(:,:,ind);
				end
				t_filterMatrices= t_filterMatrices+t_filterSecond;
			end
		end
		
		function [m_signalEstimates] = ...
				estimate(obj,noisySignal,m_basisSubspace,graph,num_localExchanges)
			%
			% Signal estimation through subspace projection.
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
			num_nodes = size(m_basisSubspace,1);
			
			% Obtain filter matrices
			[ t_filterMatrices,m_shift ] = obj.getFilterMatrices(m_basisSubspace,graph,num_localExchanges);
			
			% Signal estimation
			m_signalEstimates = NaN(num_nodes,num_localExchanges);
			for ind_localExchanges = 1:num_localExchanges
				m_filter_now = t_filterMatrices(:,:,ind_localExchanges);
				m_signalEstimates(:,ind_localExchanges) = m_filter_now*noisySignal;
			end
			
		end
		
		function [m_coefficients,normalized_error,m_filter] = findCoefficientsForOrder(obj,m_shift,m_projection,order)
			
			%
			%  m_coefficients is a size(m_shift,1) x order matrix with the
			%  node-dependent coefficients.
			%
			%  normalized_error = norm( H - P ,'fro')^2 / norm( P , 'fro' )^2
			%
			[m_coefficients,normalized_error,m_filter] = obj.findCoefficientsForOrder_static(obj.b_nodeVariant,m_shift,m_projection,order);
		
			
		end
		
		function [m_cons]=getConstraintMatrix(obj,graph)
			% it Obtains a matrix that satisfies the constraints of the problem.
			%
			% INPUT:
			%   graph : object of class Graph
			% OUTPUT:
			%   m_cons : num_nodes^2 x num_nodes^2 x
			%       a matrix that satisfies the constraints of the problem
			[b_sym]=graph.CheckSym();
			[m_connected]=graph.getConnectedEdges();
			[m_missingEdges]=graph.getMissingEdges();
			num_zero=size(m_missingEdges,2);
			num_nodes = graph.getNumberOfNodes();
			m_cons=zeros(num_nodes^2,num_nodes^2);
			if b_sym==1
				for ind_it=1:num_nodes^2-num_zero
					if m_connected(1,ind_it)~=m_connected(2,ind_it)
						m_cons(num_nodes*(m_connected(2,ind_it)-1)+m_connected(1,ind_it),num_nodes*(m_connected(2,ind_it)-1)+m_connected(1,ind_it))=.5;
						m_cons(num_nodes*(m_connected(2,ind_it)-1)+m_connected(1,ind_it),num_nodes*(m_connected(1,ind_it)-1)+m_connected(2,ind_it))=.5;
					else
						m_cons(num_nodes*(m_connected(2,ind_it)-1)+m_connected(1,ind_it),num_nodes*(m_connected(2,ind_it)-1)+m_connected(1,ind_it))=1;
					end
				end
			else
				for ind_it=1:num_nodes^2-num_zero
					m_cons(num_nodes*(m_connected(2,ind_it)-1)+m_connected(1,ind_it),num_nodes*(m_connected(2,ind_it)-1)+m_connected(1,ind_it))=1;
				end
			end
		end
		
	end
	
	methods(Static)
		
		function [U_par,U_perp,m_projection] = findOrthonormalBases(m_basisSubspace)
			% Orthonormal bases for the subspace and its orhtogonal complement
			% - document this header
			
			num_nodes = size(m_basisSubspace,1);
			dim_subspace = size(m_basisSubspace,2);
			
			U_completed=[m_basisSubspace,randn(num_nodes,num_nodes-dim_subspace)];
			[U,~]=qr(U_completed);
			U_par=U(:,1:dim_subspace);
			U_perp=U(:,dim_subspace+1:num_nodes);
			if nargout>2
				m_projection = U_par*U_par';
			end
		end
		
		function t_A = computeShiftMonomials(m_shift,max_order)
			% t_A is a size(m_shift,1) x size(m_shift,1) x max_order+1 with the
			% normalized powers of the matrix m_shift.
			
			num_nodes = size(m_shift,1);
			t_A = NaN(num_nodes,num_nodes,max_order+1);
			t_A(:,:,1) = eye(num_nodes);
			
			% Normalization (improve conditioning)
			v_evals = eig((m_shift+m_shift')/2);
			if (max(v_evals) - min(v_evals) ) > 1e-6
				m_shift = 2*(m_shift - min(v_evals)*eye(size(m_shift)))/(max(v_evals)-min(v_evals))-eye(size(m_shift));
			end
			
			
			t_A(:,:,1) = t_A(:,:,1)*(1/norm( t_A(:,:,1) , 'fro' ));
			
			for ind_order = 2:max_order+1
				t_A(:,:,ind_order) = t_A(:,:,ind_order-1)*m_shift;
				
				%normalization
				t_A(:,:,ind_order) = t_A(:,:,ind_order)*(1/norm(t_A(:,:,ind_order),'fro'));
			end
			
		end
		
		function [m_coefficients,normalized_error,m_filter] = findCoefficientsForOrder_static(b_nodeVariant,m_shift,m_target,order)
		
			
			num_nodes = size(m_shift,1);
			if length(order)>1
				error('not implemented')
			end
			if order > num_nodes - 1
				error('order must be smaller than number of nodes')
			end
			
			t_monomials = DecentralizedProjectionEstimator.computeShiftMonomials(m_shift,order);
			
			m_coefficients = NaN( num_nodes , order + 1 );
			m_filter = NaN( num_nodes, num_nodes);
			
			if b_nodeVariant
				err = zeros(1,num_nodes);
				for ind_nodes=1:num_nodes
					m_monomialsThisNode = permute(t_monomials(:,ind_nodes,:),[1 3 2]);
					v_b = m_target(:,ind_nodes);
					
					m_coefficients(ind_nodes,:) = (pinv(m_monomialsThisNode)*v_b)';
					
					err(ind_nodes) = norm(m_monomialsThisNode*m_coefficients(ind_nodes,:)'-v_b,2)^2;
					m_filter(:,ind_nodes) = m_monomialsThisNode*m_coefficients(ind_nodes,:)';
					
				end
				normalized_error_per_node = err / norm( m_target ,'fro')^2;
				normalized_error = sum(normalized_error_per_node);
			else
				m_A = reshape(t_monomials,[num_nodes^2,1,order+1]);
				m_A = permute(m_A,[1 3 2]);
				v_b = m_target(:);
				v_coefficients = pinv(m_A)*v_b;
				m_coefficients = repmat(v_coefficients,[1 num_nodes]);
				m_filter = reshape(m_A*v_coefficients,[num_nodes,num_nodes]);
				normalized_error = norm( m_filter - m_target , 'fro')^2 / norm( m_target , 'fro')^2;
				
			end
			
			% debug
			if 0 && order == num_nodes-1
				keyboard
			end
	end
		
		
	end
	
end

