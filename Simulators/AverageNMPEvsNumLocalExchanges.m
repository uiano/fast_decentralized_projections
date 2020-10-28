function [ v_error] = AverageNMPEvsNumLocalExchanges( num_MCIterations, graphGenerator , num_localExchanges , inputEstimator)
% Description: MonteCarlo simulation of the projection matrix estimation which must be
% implemented in a decentralized fashion over a topology generated using the Erdos-Renyi
% random graph model.
%
% =======>>>  UPDATE THE HEADER
% INPUT:
%   graphGenerator : object inheriting from class GraphGenerator
%   dim_subspace : dimension of the subspace containing the data
%   num_localExchanges : number of times a node shares its estimate with its
%              neighbors.
%   estimator : is an object that inherits from
%               DecentralizedProectionEstimator
%
% OUTPUT:
%   v_err : 1 x numLocalExchanges vector where the i-th entry corresponds to
%          the error between the projection matrix and the proposed filters after ind_localExchanges-th local exchanges

m_NMSENumerator = NaN(num_MCIterations,num_localExchanges);
m_NMSEDenominator = NaN(num_MCIterations,num_localExchanges);
for ind_MCIterations = 1:num_MCIterations
	
	% Generate the graph
	graph = graphGenerator.realization();
		
	% Generate the subspace and its orthogonal subspace
 	num_nodes=graph.getNumberOfNodes();
 	m_basisSubspace=(1/sqrt(num_nodes))*ones(num_nodes,1);
 	m_projection=(1/num_nodes)*ones(num_nodes);
	
	[ t_filterMatrices ] = inputEstimator.getFilterMatrices(m_basisSubspace,graph,num_localExchanges);
	
	for ind_localExchanges=1:num_localExchanges
		m_NMSENumerator(ind_MCIterations,ind_localExchanges) = norm( m_projection - t_filterMatrices(:,:,ind_localExchanges) , 'fro' )^2;
		m_NMSEDenominator(ind_MCIterations,ind_localExchanges) =  norm( m_projection , 'fro' )^2;
	end

end

v_error = mean(m_NMSENumerator,1)./mean(m_NMSEDenominator,1);
	
end