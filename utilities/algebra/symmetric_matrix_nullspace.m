function m_output = symmetric_matrix_nullspace(size_mat)
% m_output is a full rank N*(N-1)/2 x N^2 matrix with 0s, 1s, and -1s, 
% such that W*vec(A) = 0 iff A is symmetric. Can be used to enforce 
% symmetry constraints. 

num_rows = (size_mat*(size_mat-1))/2;
assert(size_mat>=2);

% Obtain pairs of indices (i,j) with i<j between 1 and size_mat
m_i = (1:size_mat)'*ones(1,size_mat);
m_j = m_i';
v_sym=triu(ones(size_mat),1);
v_upper_triang_indicator = v_sym(:)>0;
m_ij_pairs = [m_i(v_upper_triang_indicator) m_j(v_upper_triang_indicator) ];

%Obtain indices of 1's
m_indices_p1 = [ (1:num_rows)' , m_ij_pairs(:,1)+size_mat*(m_ij_pairs(:,2)-1) ];
m_indices_m1 = [ (1:num_rows)' , m_ij_pairs(:,2)+size_mat*(m_ij_pairs(:,1)-1) ];

% Output matrix
m_output = zeros( num_rows , size_mat^2 );
m_output( m_indices_p1(:,1) + num_rows*(m_indices_p1(:,2)-1) )  = 1;
m_output( m_indices_m1(:,1) + num_rows*(m_indices_m1(:,2)-1) )  = -1;

end 
