function [ v_error] = FilterNMSEvsNumLocalExchanges( num_MCIterations, graphGenerator , dim_subspace , num_localExchanges , inputEstimator, num_failures)
% Description: MonteCarlo simulation of the projection matrix estimation which must be
% implemented in a decentralized fashion over a topology generated using the Erdos-Renyi
% random graph model.
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
if nargin < 6
	num_failures = 0
end

m_NMSENumerator = NaN(num_MCIterations,num_localExchanges);
m_NMSEDenominator = NaN(num_MCIterations,num_localExchanges);
for ind_MCIterations = 1:num_MCIterations
	
	% Generate the graph
	graph = graphGenerator.realization();
	% Generate the subspace and its orthogonal subspace
	[m_basisSubspace,~]=qr(rand(graph.getNumberOfNodes(),dim_subspace),0);
	m_projection = m_basisSubspace*m_basisSubspace';
	
	if num_failures
		graph = graph.applyNodeFailures(num_failures);
	end
	
	[ t_filterMatrices] = inputEstimator.getFilterMatrices(m_basisSubspace,graph,num_localExchanges);
	
	for ind_localExchanges=1:num_localExchanges
		m_NMSENumerator(ind_MCIterations,ind_localExchanges) = norm( m_projection - t_filterMatrices(:,:,ind_localExchanges) , 'fro' )^2;
		m_NMSEDenominator(ind_MCIterations,ind_localExchanges) =  norm( m_projection , 'fro' )^2;
	end
	
end

v_error = mean(m_NMSENumerator,1)./mean(m_NMSEDenominator,1);

end
