function [ prob_feasible] = ProbFeasibilityvsNumLocalExchanges( num_MCIterations, graphGenerator , dim_subspace , estimator)
% Description: This function uses Monte Carlo simulation to estimate the probability that
% the problem of obtaining a valid shift matrix is feasible. 
%
% INPUT:
%   graphGenerator : object inheriting from class GraphGenerator
%   dim_subspace : dimension of the subspace containing the data
%   estimator : is an object that inherits from
%               DecentralizedProectionEstimator. It must implement the
%               method isFeasible().
%
% OUTPUT:
%   prob_feasible : estimate of the probability that the problem is
%   infeasible

v_feasibleHistory = NaN(1,num_MCIterations);
for ind_MCIterations = 1:num_MCIterations
	
	% Generate the graph
	graph = graphGenerator.realization();
		
	% Generate the subspace and its orthogonal subspace
	[m_basisSubspace,~]=qr(rand(graph.getNumberOfNodes(),dim_subspace),0);
		
	v_feasibleHistory(ind_MCIterations) = estimator.isFeasible(m_basisSubspace,graph);
	
end

prob_feasible = mean(v_feasibleHistory);
	
end
