classdef DLMSEDecentralizedProjectionEstimator < DecentralizedProjectionEstimator
	
	properties
		
		stepSizeDLMS %the step size used in the x-minimization iteration of the ADMM method
		lagrangeDLMS %the augmented Lagrangian parameter
	end
	
	methods
		function [m_signalEstimatesDLMS]=estimate(obj,noisySignal,m_basisSubspace,graph,num_localExchanges)
			
			% input processing
			num_nodes = size(m_basisSubspace,1);
			subspaceDimension=size(m_basisSubspace,2);
			
			% Initialization matrices
			m_lagrange=zeros(subspaceDimension,num_nodes,num_nodes);% Matrix of the lagrange maltipliers of ADMM
			m_object=zeros(subspaceDimension,num_nodes);% the signal that we want to estimate (the object of the estimation)
			
			% obtaining the connected edges
			m_connected=graph.getConnectedEdges();
			
			% program input
			m_signalEstimatesDLMS=zeros(num_nodes,num_localExchanges);
			
			% program
			for ind_numDLMSIterate=1:num_localExchanges
				for ind_numnodes=1:num_nodes
					
					%input initialization
					E_1=zeros(subspaceDimension,1);
					E_2=zeros(subspaceDimension,1);
					
					%ADMM updates
					[rowLocation,columnLocation]=find(m_connected==ind_numnodes); %finding the nodes which are connected to ind_numnodes-th node by using m_connected
					for lengthLocation=1:size(columnLocation,1)
						if rowLocation(lengthLocation)==1
							E_1=E_1+m_lagrange(:,m_connected(columnLocation(lengthLocation)),ind_numnodes)-...%it is needed for the updates of ADMM
								m_lagrange(:,ind_numnodes,m_connected(columnLocation(lengthLocation)));
							E_2=E_2+m_object(:,ind_numnodes)-(m_object(:,m_connected(columnLocation(lengthLocation))));%it is needed for the updates of ADMM
						end
					end
					m_object(:,ind_numnodes)=m_object(:,ind_numnodes)+(obj.stepSizeDLMS)*(2*m_basisSubspace(ind_numnodes,:)'...%update of the object signal based on the ADMM updates
						*(noisySignal(ind_numnodes)-m_basisSubspace(ind_numnodes,:)*m_object(:,ind_numnodes))-E_1-obj.lagrangeDLMS*E_2);
				end
				
				for ind_numnodes=1:num_nodes
					[rowLocation,columnLocation]=find(m_connected==ind_numnodes);
					for lengthLocation=1:size(columnLocation,1)
						if rowLocation(lengthLocation)==1
							m_lagrange(:,m_connected(columnLocation(lengthLocation)),ind_numnodes)=...% update of the lagrange multipliers based on the ADMM updates
								m_lagrange(:,m_connected(columnLocation(lengthLocation)),ind_numnodes)...
								+(obj.lagrangeDLMS/2)*(m_object(:,ind_numnodes)-(m_object(:,m_connected(columnLocation(lengthLocation)))));
						end
					end
				end
				v_signalEstimatesDLMS=NaN(1,num_nodes);
				for ind_nodes=1:num_nodes
					v_signalEstimatesDLMS(ind_nodes)=m_basisSubspace(ind_nodes,:)*m_object(:,ind_nodes);
				end
				m_signalEstimatesDLMS(:,ind_numDLMSIterate)= v_signalEstimatesDLMS';%the signal estimation
			end
		end
	end
end