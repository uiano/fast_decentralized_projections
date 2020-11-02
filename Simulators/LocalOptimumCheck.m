function [ v_error] = LocalOptimumCheck( graphGenerator , dim_subspace , num_localExchanges ,v_alpha, inputEstimator,v_gamma,b_control)
% Description: MonteCarlo simulation to check whether the Edge-variant and
%              the proposed method find a local optimum or not.
% INPUT:
%   num_MCIterations: The number of Montecarlo simulation
%   graphGenerator : object inheriting from class GraphGenerator
%   dim_subspace : dimension of the subspace containing the data
%   num_localExchanges : number of times a node shares its estimate with its
%              neighbors.
%    v_alpha: The vector contains alpha
%   estimator : is an object that inherits from
%               DecentralizedProectionEstimator
%     b_control: the proposed method checks if  b_control=1
%                the Edge-Variant method is checked if v_control=0
%     v_gamma:  the vector contains regularization coefficients of
%                the succesive method  
% OUTPUT:
%   v_err : 1 x numLocalExchanges vector where the i-th entry corresponds to
%          the error between the projection matrix and the proposed filters after ind_localExchanges-th local exchanges
% Generate the graph
graph = graphGenerator.realization();
m_adjacency=graph.m_adjacency;
% Generate the subspace and its orthogonal subspace
[m_basisSubspace,~]=qr(rand(graph.getNumberOfNodes(),dim_subspace),0);
m_projection = m_basisSubspace*m_basisSubspace';
num_nodes=graph.getNumberOfNodes();
[ t_filterMatrices,m_shift] = inputEstimator.getFilterMatrices(m_basisSubspace,graph,num_localExchanges);
for ind_num=1:length(v_alpha)
	for ind=1:num_localExchanges
		m_rand(:,:,ind)=rand(num_nodes,num_nodes);
		S_new(:,:,ind)=m_shift(:,:,ind)+v_alpha(ind_num)*m_rand(:,:,ind);
		S_new(:,:,ind)=S_new(:,:,ind).*(m_adjacency>0);
	end
	m_filter=eye(num_nodes,num_nodes);
	m_agg=zeros(num_nodes,num_nodes,num_localExchanges);
	m_sum=zeros(num_nodes,num_nodes);
	if b_control==0
		for i=1:num_localExchanges
			m_filter=S_new(:,:,i)*m_filter;
			m_sum=m_sum+m_filter;
			m_agg(:,:,i)=m_sum;
		end
	else
		for i=1:num_localExchanges
			m_filter=S_new(:,:,i)*m_filter;
			m_agg(:,:,i)=m_filter;
		end
	end
	m_sum2=0;
	if b_control==0
		v_error(ind_num)=norm(m_projection-m_agg(:,:,end),'fro')^2;
	else
		for ind=1:length(v_gamma)
			m_sum2=m_sum2+v_gamma(ind)*norm(m_projection-m_agg(:,:,ind),'fro')^2;
		end
		v_error(ind_num)=m_sum2;
	end
end