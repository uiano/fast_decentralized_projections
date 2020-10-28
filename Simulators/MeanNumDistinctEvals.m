function [ num_evals,v_eigPar,v_eigPerp] = MeanNumDistinctEvals( num_MCIterations , graphGenerator , dim_subspace, evalTol, inputEstimator)
% Description: Monte Carlo simulation for obtaining the number of
% distinct eigenvalues.
%
%
% INPUT:
%   num_MCiterations : number of Monte Carlo iterations
%   graphGenerator  : an object of class Graph
%   dim_subspace : dimension of the subspace containing the data
%   evalTol :  two eigenvalues are equal whenever the absolute value of
%    their difference is less than or equal to evalTol.
%   estimator : is an obect that inherits from
%               DecentralizedProectionEstimator
%
% OUTPUT:
%   num_evals : The number of distinct eigenvalues

%input processing
num_lambdaPar=NaN(1,num_MCIterations);
num_lambdaPerp=NaN(1,num_MCIterations);
num_differEig=NaN(1,num_MCIterations);

% Monte Carlo loop
for ind_MCIteration = 1:num_MCIterations
	% Generate the graph
	graph = graphGenerator.realization();
	
	% Generate the subspace and its orthogonal subspace
	[m_basisSubspace,~]=qr(rand(graph.getNumberOfNodes(),dim_subspace),0);
	
	%obtaining F_par (in the code is m_par) and F_perp (in the code is m_perp)
	[~,~,m_shiftPar,m_shiftPerp] = inputEstimator.findShiftMatrix(m_basisSubspace,graph);
	[U_par,U_perp] = inputEstimator.findOrthonormalBases(m_basisSubspace);
	m_par=U_par'*m_shiftPar*U_par;
	m_perp=U_perp'*m_shiftPerp*U_perp;
	
	
	%Obtaining the number of distinct eigenvalues
	[num_lambdaPar(ind_MCIteration),v_eigPar]=DistinctEvals(m_par,evalTol);
	[num_lambdaPerp(ind_MCIteration),v_eigPerp]=DistinctEvals(m_perp,evalTol);
	num_differEig(ind_MCIteration)=num_lambdaPerp(ind_MCIteration)+num_lambdaPar(ind_MCIteration);
end
% Obtain the average of the number of distinct eigenvalues
num_evals=mean(num_differEig);

end
function  [num_lambda,v_eig]=DistinctEvals(m_input,evalTol)
v_eig=eig(m_input)';%eigenvalues of m_par
v_lambda=UniqueTolerance(v_eig,evalTol);%finding the unique eigenvalues of m_input
num_lambda=length(v_lambda);
end
function [lambda_unique]=UniqueTolerance(lambda,evalTol)
% Description: It finds the unique eigenvalues

while 1
	lambda=sort(lambda);
	lambda_differe=lambda(2:length(lambda))-lambda(1:length(lambda)-1);
	[differe_min, position]=min(lambda_differe);
	if abs(differe_min)<evalTol
		lambda(position)=(lambda(position)+lambda(position+1))/2;
		lambda(position+1)=[];
	else
		break
	end
end
lambda_unique=unique(lambda);
end