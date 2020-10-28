function [m_adjacency]=graphGenerator(N,missingEdgeProbability)
m_adjacency = ( rand(N,N)>= missingEdgeProbability );
m_adjacency= m_adjacency - tril(m_adjacency,-1) + triu(m_adjacency,1)';
m_adjacency(logical(eye(N))) = 1;
end