classdef DGDDecentralizedProjectionEstimator < DecentralizedProjectionEstimator
	
	properties
		
		stepSizeDGD %the step size used for the gradient method
	end
	
	methods
		function [m_signalEstimatesDGD]=estimate(obj,noisySignal,m_basisSubspace,graph,num_localExchanges)
			
			% input processing
			num_nodes = size(m_basisSubspace,1);
			subspaceDimension=size(m_basisSubspace,2);
			
			
			% Initialization matrices
			m_coordinates=zeros(subspaceDimension,num_nodes);% each column of m_coordinates contains the estimation of the signalCoordinates that belongs to each node
			m_counter=NaN(subspaceDimension,num_nodes);%saves sum_j=1^num_nodes m_adjacencyDGD_ij*x_j^k at its columns where k is the number of the iterations
			% and x is a signal that the DGD method wants to estimate it
			m_signalEstimatesDGD=zeros(num_nodes,num_localExchanges);%the goal of estimation m_signalEstimatesDGD=m_basisSubspace*m_coordinates. Indeed, we have
			%noisysignal=m_basisSubspace*signalCoordinates+noise. The DGD method estimates signalCoordinates called as m_coordinates. Then we can find the estimation of
			%m_basisSubspace*signalCoordinates called as m_signalEstimatesDGD which we want.
			
			%obtaining the adjacency matrix for DGD
			m_adjacencyDGD=graph.m_adjacency./sum(graph.m_adjacency,2);%this matrix guarantees that we will reach a consensus among nodes finally. Indeed, in
			% the constraint of a problem which the DGD method solvs we have: m_adjacencyDGD*x=x
			
			%program
			for ind_numDGDIterate=1:num_localExchanges
				for ind_node=1:num_nodes
					m_counter(:,ind_node)=(m_adjacencyDGD(ind_node,:)*m_coordinates')';%saves sum_j=1^num_nodes m_adjacencyDGD_ij*x_j^k where k is the number of the iterations
					% and x is a signal that we want to estimate at m_counter columns
					m_coordinates(:,ind_node)=m_counter(:,ind_node)-((obj.stepSizeDGD))*(2*m_basisSubspace(ind_node,:)'*...
						m_basisSubspace(ind_node,:)*m_coordinates(:,ind_node)-2*m_basisSubspace(ind_node,:)'*noisySignal(ind_node));...
						%updates the estimated signal (x_i) at each node i.e x_i^k+1=sum_j=1^num_nodes W_ij*x_j^k-stepsize*gradient of frobeniuse norm(x_i^k)
				end
				v_signalEstimatesDGD=NaN(1,num_nodes);
				for ind_nodes=1:num_nodes
					v_signalEstimatesDGD(ind_nodes)=m_basisSubspace(ind_nodes,:)*m_coordinates(:,ind_nodes);
				end
				m_signalEstimatesDGD(:,ind_numDGDIterate)= v_signalEstimatesDGD';%the estimated signal 
			end
		end
	end
end