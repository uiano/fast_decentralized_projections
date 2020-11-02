classdef FastDecentralizedProjectionEstimator < DecentralizedProjectionEstimator
	% Description: This class represents our proposed method which enables
	% distributed subspace projection via graph filters.
	
	properties
		
		%% GENERAL PROPERTIES
		
		str_criterion = 'exact';
		%   Criterion to find the optimum shift. Possible values:
		%       'exact' : Attempts to find a graph filter to compute the
		%                 projection exactly while minimizing the number of
		%                 local exchanges
		%       'approximate' : attempts to find a graph filter that
		%                 approximates the projection and minimizes the
		%                 number of local exchanges
		%        'auto' : chooses automatically between 'exact' and
		%                 'approximate'
		%
		epsilon = NaN;
		%     This property pertains to the constraint
		%     tr(F_perp)=(N-r)(1-epsilon) which prevents the trivial solutions
		%     If epsilon = NaN, then such a constraint is not enforced.
		verboseLevel = 0;
		%     Set this variable to a positive integer to obtain additional information during
		%     run time. The higher the integer, the more information.
		maxParameterValues = 10;
		%     Maximum number of parameter values to try when tuning.
		b_normalizeRegParameters = 1; % If 1, the regularization parameters are divided by their sum.
		b_nodeVariant = 1; % set to 0 to use node-invariant filters
		str_method = 'ADMM';
		num_alt=10;
		num_layer;
		%This property denotes that which method should work (ADMM or PGD)
		%% REGULARIZATION PARAMETERS
		%  When str_criterion = 'exact':
		%     obj.reg_evalMultiplicityPar * nuc_kron( F_par )  / dim_subspace^2
		%     + obj.reg_evalMultiplicityPerp * nuc_kron( F_perp ) / (num_nodes - dim_subspace)^2
		%     + obj.reg_energyPar * frobenius_square( F_par )  / dim_subspace^2
		%     + obj.reg_energyPerp * frobenius_square( F_perp ) / (num_nodes - dim_subspace)^2
		% When str_criterion = 'approximate':
		%     obj.reg_evalMultiplicityPar * nuc_kron( U_par' * S * U_par ) / dim_subspace^2
		%     + obj.reg_evalMultiplicityPerp * nuc_kron( U_perp' * S * U_perp  )  / (num_nodes - dim_subspace)^2
		%     + obj.reg_separation * frobenius_square( U_par' * S * U_perp )  / ( (num_nodes - dim_subspace) dim_subspace )
		%     + obj.reg_energyPar * frobenius_square( U_par' * S * U_par  ) / dim_subspace^2
		%     + obj.reg_energyPerp * frobenius_square( U_perp' * S * U_perp )  / (num_nodes - dim_subspace)^2
		reg_evalMultiplicityPar% = 1;
		reg_evalMultiplicityPerp% = 0;
		reg_separation% = 0.01;
		reg_energyPar% = 1;
		reg_energyPerp% = 0;
		toleranceFilterError = .005; % the function autoSetRegParameters()
		%  will attempt to find the combination regularization parameters
		%  for which an error norm( m_filter - m_projection , 'fro')^2 /
		%  num_nodes^2 less than toleranceFilterError is attained with
		%  (ideally) the smallest number of local exchanges.
		
		%% PROPERTIES FOR THE EXACT PROJECTION CRITERION
		num_maxIterationsPSGD = 1000;
		stepSizePSGD = 50; %% PSGD uses stepSizePSGD/ind_iteration as step size
		tolerancePSGD = 1e-7; %% if
		%  norm( previous iterate - new iterate )/( dim_subspace^2 + (num_nodes-dim_subspace)^2 ) < tolerancePSGD,
		%  then PSGD stops.
		
		%% PROPERTIES FOR THE APPROXIMATE PROJECTION CRITERION
		num_maxIterationsADMM=1000
		% the number of ADMM iterations
		stepSizeADMM = .1;
		%reg_firstAugmented=1
		%reg_secondAugmented=0.5
		toleranceADMM = 1e-7;%% if norm( previous iterate - new iterate , 'fro' )/num_nodes^2 < toleranceADMM, ADMM stops.
	     v_gamma; % the weighting vector
	end
	
	
	methods(Access = public)
		
		function [m_S,m_projection,m_shiftPar,m_shiftPerp,v_objective] = findShiftMatrix(obj,m_basisSubspace,graph)
			
			% Parameter check
			assert(obj.reg_evalMultiplicityPar >= 0);
			assert(obj.reg_evalMultiplicityPerp >= 0);
			assert(obj.reg_energyPar >= 0);
			assert(obj.reg_energyPerp >= 0);
			
			% Orthonormal bases
			[U_par,U_perp,m_projection] = obj.findOrthonormalBases(m_basisSubspace);
			
			
			% Obtain shift matrix
			switch obj.str_method
				
				%obtain the shift matrix of PGD
				case('PGD')
					[m_S,m_shiftPar,m_shiftPerp,~,~,v_objective] = obj.exactShiftSolver(U_par,U_perp,graph,obj.epsilon);
					
				%Obtain the shift matrix of CVX
				case('CVX')
						[m_S]=obj.CVXShiftSolver(U_par,U_perp,graph);
					
				case('ADMM')
					switch	obj.str_criterion
						case('exact')							
							
							% Reg. parameter normalization
							if obj.b_normalizeRegParameters
								regParSum = obj.reg_evalMultiplicityPar + obj.reg_evalMultiplicityPerp ...
									+ obj.reg_energyPar  + obj.reg_energyPerp ;
								obj.reg_evalMultiplicityPar = obj.reg_evalMultiplicityPar / regParSum;
								obj.reg_evalMultiplicityPerp  = obj.reg_evalMultiplicityPerp / regParSum;
								obj.reg_separation = obj.reg_separation / regParSum;
								obj.reg_energyPar = obj.reg_energyPar /regParSum;
								obj.reg_energyPerp = obj.reg_energyPerp / regParSum;
							end
							
							if (obj.reg_evalMultiplicityPar == 0) && (obj.reg_evalMultiplicityPerp == 0 )
								% Objective does not include nuclear norms
								[m_S,m_shiftPar,m_shiftPerp] =	obj.exactShiftSolverNoFast(U_par,U_perp,graph,obj.epsilon);
							else
								% Objective includes nuclear norms
								[m_S,m_shiftPar,m_shiftPerp,~,~,v_objective] = obj.exactShiftSolverADMM(U_par,U_perp,graph,obj.epsilon);
							end
							
						case('approximate')
							
							assert(obj.reg_separation >= 0);
							
							% Reg. parameter normalization
							if obj.b_normalizeRegParameters
								regParSum = obj.reg_evalMultiplicityPar + obj.reg_evalMultiplicityPerp ...
									+ obj.reg_separation + obj.reg_energyPar  + obj.reg_energyPerp ;
								obj.reg_evalMultiplicityPar = obj.reg_evalMultiplicityPar / regParSum;
								obj.reg_evalMultiplicityPerp  = obj.reg_evalMultiplicityPerp / regParSum;
								obj.reg_separation = obj.reg_separation / regParSum;
								obj.reg_energyPar = obj.reg_energyPar /regParSum;
								obj.reg_energyPerp = obj.reg_energyPerp / regParSum;
							end
							
							% Shift computation
							if (obj.reg_evalMultiplicityPar == 0) && (obj.reg_evalMultiplicityPerp == 0 )
								% Objective does not include nuclear norms
								[m_S]=	obj.approximateShiftSolverNoFast(U_par,U_perp,graph);
							else
								% Objective includes nuclear norms
								[m_S]= obj.approximateShiftSolver(U_par,U_perp,graph,obj.epsilon);
							end
							
						case('auto')
							error('not implemented')
							
					end
			end
			
			m_S = m_S.*(graph.m_adjacency~=0);
			
		end
		
		function [t_filterMatrices,m_shift]=findSecondShift(obj,m_differ,graph)
		    % This method belongs to the successive paper
			% it Obtains a filter and the graph shift operator by using the successive method).
			% Indeed, in this method two GFs is used (Fast+Successive)
			% INPUT:
			%  m_differ: The difference between the first filter
			%  and m_projection
			%   graph : object of class Graph
			% OUTPUT:
			%   t_filterMatrices : num_nodes x num_nodes x
			%       num_local_exchanges tensor where
			%       t_filterMatrices(:,:,n) is the filter matrix
			%       after n local exchanges.
			%   m_shift : num_nodes x num_nodes x
			%       matrix where m_shift
			%        is the shift matrix
			num_nodes=size(m_differ,1);
			m_shift=zeros(num_nodes,num_nodes,obj.num_layer);
			for ind_s=1:obj.num_layer
				m_shift(:,:,ind_s)=eye(num_nodes,num_nodes);
			end
			[m_cons]=obj.getConstraintMatrix(graph);
			% Algorithm
			for i=1:obj.num_alt
				for ind_deep=1:obj.num_layer
					counter=0;
					m_mulleft=eye(num_nodes,num_nodes);
					m_mulright=eye(num_nodes,num_nodes);
					m_counter=NaN(num_nodes,num_nodes,obj.num_layer-ind_deep+1);
					m_counter(:,:,1)=eye(num_nodes);
					if ind_deep<obj.num_layer
						for ind_mul=ind_deep+1:1:obj.num_layer
							m_mulleft=m_shift(:,:,ind_mul)*m_mulleft;
							counter=counter+1;
							m_counter(:,:,counter+1)=m_mulleft;
						end
					end
					if ind_deep>1
						for ind_mul2=ind_deep-1:-1:1
							m_mulright=m_mulright*m_shift(:,:,ind_mul2);
						end
					end
					m_inv=zeros(num_nodes^2,num_nodes^2);
					m_sec=zeros(num_nodes^2,num_nodes^2);
						count=0;
						for ind_cost=ind_deep:obj.num_layer
							count=count+1;
							m_frist=kron(m_mulright',(m_counter(:,:,count)))*m_cons;
							m_inv=obj.v_gamma(ind_cost)*(m_frist)'*m_frist+m_inv;
							m_sec=obj.v_gamma(ind_cost)*m_frist'+m_sec;
						end
					S=pinv(m_inv)*(m_sec*m_differ(:));
					m_shift(:,:,ind_deep)=reshape(S,num_nodes,num_nodes);
					m_shift(:,:,ind_deep) =m_shift(:,:,ind_deep).*(graph.m_adjacency>0);
				end
			end
			m_filter=eye(num_nodes,num_nodes);
			m_agg=zeros(num_nodes,num_nodes,obj.num_layer);
			for i=1:obj.num_layer
				m_filter=m_shift(:,:,i)*m_filter;
				m_agg(:,:,i)=m_filter;
			end
			t_filterMatrices=m_agg;
		end
		
		function b_feasible = isFeasible(obj,m_basisSubspace,graph,epsilon_in)
			%
			% Returns 1 if the optimization solved to find the shift matrix
			% is feasible; 0 otherwise.
			%
			
			% input processing
			if nargin < 4
				epsilon_in = obj.epsilon;
			end
			
			[U_par,U_perp] = obj.findOrthonormalBases(m_basisSubspace);
			
			[~,~,b_feasible] = obj.feasibleSetEquations(U_par,U_perp,graph,epsilon_in);
		end
		
		
		function [m_A,v_b,m_size,b_isFeasible,dim_prefeasible] = feasibleSetEquations(obj,U_par,U_perp,graph,epsilon)
			%
			% Provides a matrix m_A and vector v_b such that the set of
			% feasible shift matrices is given by the solution to the
			% linear system of equations
			%
			%        m_A * v_x  = v_b
			%
			% where, depending on obj.str_criterion:
			%
			%  - if obj.str_criterion == 'exact': then the solutions
			%    v_x = [ F_par(:);F_perp(:)] correspond to the feasible
			%    shifts for the exact projection problem.
			%
			%  - if obj.str_criterion == 'approximate': then the solutions
			%    v_x = [ m_shift(:) ] correspond to the feasible
			%    shifts for the approximate projection problem.
			%
			% In both cases, the trace constraint with (1+obj.epsilon)
			% is used.
			%
			% INPUT:
			%   U_par: num_nodes x dim_subspace matrix
			%       whose orthonormal columns span the signal subspace.
			%   U_perp: num_nodes x (num_nodes - dim_subspace) matrix
			%       whose orthonormal columns span the orthogonal complement
			%       to the signal subspace.
			%   graph : object of class Graph
			%   epsilon [optional]: value of epsilon in the trace
			%       constraint; default is obj.epsilon.
			% OUTPUT:
			%   b_isFeasible: 1 if feasible set not empty; 0 otherwise.
			%   dim_feasibleSet : dimension of the subspace of prefeasible
			%       shift matrices.
			%   m_A and v_b : explained above.
			%
			
			% input processing
			num_nodes = size(U_par,1);
			dim_subspace = size(U_par,2);
			if nargin<5
				epsilon = obj.epsilon;
			end
			
			% Find matrices common to both exact and approx problem
			[m_missingEdges] = graph.getMissingEdges();
			[m_topologyNullSpace] = obj.symmetricShiftNullSpace(num_nodes,m_missingEdges);
			m_size=size(m_topologyNullSpace);
			
			switch obj.str_criterion
				case('exact')
					
					% form a matrix whose nullspace determines the set of prefeasible
					% shift matrices
					if dim_subspace>1
					m_symmetryF_par = symmetric_matrix_nullspace(dim_subspace);
					m_symmetryF_perp = symmetric_matrix_nullspace(num_nodes - dim_subspace);
					m_preFeasibleA = [m_topologyNullSpace*kron(U_par,U_par),m_topologyNullSpace*kron(U_perp,U_perp);
						m_symmetryF_par, zeros(size(m_symmetryF_par,1),(num_nodes-dim_subspace)^2 ) ;
						zeros(size(m_symmetryF_perp,1),dim_subspace^2),m_symmetryF_perp];
					else
					m_symmetryF_perp = symmetric_matrix_nullspace(num_nodes - dim_subspace);
					m_preFeasibleA = [m_topologyNullSpace*kron(U_par,U_par),m_topologyNullSpace*kron(U_perp,U_perp);
						zeros(size(m_symmetryF_perp,1),dim_subspace^2),m_symmetryF_perp];
					end
					
					% form a matrix and vector to avoid ambiguity in the solution and avoid
					% trivial solution = identity
					if dim_subspace>1
						m_trace_par = eye(dim_subspace);
						m_trace_perp = eye(num_nodes - dim_subspace);
						m_ambiguity = [ m_trace_par(:)', zeros(1,(num_nodes-dim_subspace)^2) ];
						v_ambiguity = dim_subspace;
						if size(m_missingEdges,2)/num_nodes^2>=.83
							epsilon=NaN;
						end
						if ~isnan(epsilon)
							m_ambiguity = [m_ambiguity ; zeros(1,dim_subspace^2) , m_trace_perp(:)'];
							v_ambiguity = [v_ambiguity;(1+epsilon)*(num_nodes-dim_subspace)];
						end
					else
						m_trace_perp = eye(num_nodes - dim_subspace);
						if size(m_missingEdges,2)/num_nodes^2>=.83
							epsilon=NaN;
						end
						if ~isnan(epsilon)
							m_ambiguity = [ zeros(1,dim_subspace^2) , m_trace_perp(:)'];
							v_ambiguity = [(1+epsilon)*(num_nodes-dim_subspace)];
						end
					end
					
				case('approximate')
					% form a matrix whose nullspace determines the set of prefeasible
					% shift matrices
					m_preFeasibleA = [m_topologyNullSpace;
						symmetric_matrix_nullspace(num_nodes)];
					
					% form a matrix and vector to avoid ambiguity in the solution and avoid
					% trivial solution = identity
					m_trace_par = U_par*U_par';
					m_trace_perp = U_perp*U_perp';
					m_ambiguity = [ m_trace_par(:)'];
					v_ambiguity = [dim_subspace];
					if ~isnan(epsilon)
						v_ambiguity = [v_ambiguity ;(1+epsilon)*(num_nodes-dim_subspace)];
						m_ambiguity = [m_ambiguity ; m_trace_perp(:)'];
					end
					
				otherwise
					error('obj.str_criterion must be either ''exact'' or ''approximate''');
			end
			
			% output computation
			m_A = [m_preFeasibleA;m_ambiguity];
			v_b = [zeros(size(m_preFeasibleA,1),1);v_ambiguity];
			
			if nargout > 2
				b_isFeasible = ( rank(m_A) == rank( [m_A v_b] ) );
				dim_prefeasible = size(m_preFeasibleA,2)-rank(m_preFeasibleA);
			end
			
		end
		
	end
	
	
	methods(Access = private)
		
		function [m_shift,m_perp,m_par,m_shiftPar,m_shiftPerp,v_objectivePGD] = exactShiftSolver(obj,U_par,U_perp,graph,epsilon)
			
			
			% input processing
			dim_subspace=size(U_par,2);
			num_nodes = size(U_par,1);
			
			% constants
			num_iterations = obj.num_maxIterationsPSGD;
			
			% obtaining the missing edges
			[m_missingEdges] = graph.getMissingEdges();
			
			% initialization
			m_shiftPar=randn(dim_subspace,dim_subspace);
			m_shiftPar=(m_shiftPar+m_shiftPar')/2;
			m_shiftPerp=randn(num_nodes-dim_subspace,num_nodes-dim_subspace);
			m_shiftPerp=(m_shiftPerp+m_shiftPerp')/2;
			m_shiftParOld = - m_shiftPar;
			m_shiftPerpOld = - m_shiftPerp;
			for iterationIndex=1:num_iterations
				
				if obj.verboseLevel > 5
					rp = obj.reg_evalMultiplicityPar
					iterationIndex
				end
				
				% Subgradient computation
				[m_par_nuclear]=obj.subgradientNuclearKronecker(m_shiftPar);%subgradient of nuclear_kroncker(m_shiftPar) operator
				[m_perp_nuclear]=obj.subgradientNuclearKronecker(m_shiftPerp);%subgradient of nuclear_kroncker(m_shiftPerp) operator
				m_par_subgradient=obj.reg_evalMultiplicityPar*m_par_nuclear/dim_subspace^2+obj.reg_energyPar*(2*m_shiftPar)/dim_subspace^2;%subgradient of nuclear_kron+Frobenius norm
				m_perp_subgradient=obj.reg_evalMultiplicityPerp*m_perp_nuclear/(num_nodes-dim_subspace)^2+obj.reg_energyPerp*(2*m_shiftPerp)/(num_nodes-dim_subspace)^2;%subgradient of nuclear_kron+Frobenius norm
				
				% Subgradient step
				m_par_step = m_shiftPar - obj.stepSizePSGD/iterationIndex*m_par_subgradient;
				m_per_step = m_shiftPerp - obj.stepSizePSGD/iterationIndex*m_perp_subgradient;
				
				% Projection
				[m_shiftPar,m_shiftPerp]=FastDecentralizedProjectionEstimator.projection_FF(m_par_step,U_par,m_per_step,epsilon,m_missingEdges,U_perp);%projection
				
				% Objective evaluation
				%if obj.verboseLevel
					v_objectivePGD(1,iterationIndex) = obj.objectiveExactProjection(m_shiftPar,m_shiftPerp); % speed is not important in verbose mode
				%end
				
				% Stopping criterion
				if norm( [m_shiftPar(:);m_shiftPerp(:)] - [m_shiftParOld(:);m_shiftPerpOld(:)] )/( dim_subspace^2 + (num_nodes-dim_subspace)^2 ) < obj.tolerancePSGD
					break
				end
				m_shiftParOld = m_shiftPar;
				m_shiftPerpOld = m_shiftPerp;
			end
			
			% Issue warning if not converged
			if iterationIndex == num_iterations
				warning('PSGD has not converged within the target tolerance')
				obj
			end
			
			m_par=U_par*m_shiftPar*U_par';%obtain S_parallel
			m_perp=U_perp*m_shiftPerp*U_perp';
			m_shift=m_par+m_perp;
			
			if obj.verboseLevel == 5 % compare with quadratic solver
				obj
				[~ ,F_par,F_perp]= exactShiftSolverNoFast(obj,U_par,U_perp,graph,epsilon);
				objQuadraticSolver = obj.objectiveExactProjection(F_par,F_perp)
				objPSGD = obj.objectiveExactProjection(m_shiftPar,m_shiftPerp)
				
			end
			
			% Plot objective
			if obj.verboseLevel >= 2 % extra verbose level
				
				F = GFigure('m_Y',v_objectivePGD,'ch_title','Projected Subgradient Descent',...
					'ch_ylabel','Objective value','ch_xlabel','Iteration index');
				F.plot(3);
				pause()
				
			end
			
			
			
		end
		
		
		function [m_shift,m_shiftPar,m_shiftPerp,v_perp,v_par,v_objectiveADMM] = exactShiftSolverADMM(obj,U_par,U_perp,graph,epsilon)
			%
			% m_shift is the function that minimizes the exact
			% projection criterion described in the properties section.
			%
			% The algorithm implemented here is based on the scaled form of
			% ADMM.
			%
			% This is an initial implementation. Some performance improvements are
			% still possible.
			
			% Input processing
			num_nodes = size(U_par,1);
			dim_subspace = size(U_par,2);
			rho = obj.stepSizeADMM;
			
			% Matrix for m_shift update
			m_kroneckerDifOpPar = obj.operatorMatrixKoneckerDifference( dim_subspace );
			m_kroneckerDifOpPerp = obj.operatorMatrixKoneckerDifference( num_nodes - dim_subspace );
			[m_linearCons,v_constraint] = obj.feasibleSetEquations(U_par,U_perp,graph,epsilon);
			m_linearConspar=m_linearCons(:,1:dim_subspace^2);
			m_linearConsper=m_linearCons(:,dim_subspace^2+1:end);
			
			%initialization
			v_par=randn(dim_subspace^2,1);
			v_perp=randn((num_nodes-dim_subspace)^2,1);
			v_lagrangefirst=zeros(size(v_constraint));
			v_lagrangesecond=zeros(dim_subspace^4,1);
			v_lagrangethird=zeros((num_nodes-dim_subspace)^4,1);
			m_invpar=inv(2*obj.reg_energyPar*eye(dim_subspace^2)+rho*m_linearConspar'*m_linearConspar+(rho)*m_kroneckerDifOpPar'*m_kroneckerDifOpPar);
			m_invper=inv(2*obj.reg_energyPerp*eye((num_nodes-dim_subspace)^2)+rho*m_linearConsper'*m_linearConsper+(rho)*m_kroneckerDifOpPerp'*m_kroneckerDifOpPerp);
			m_shiftParOld=- eye(num_nodes);
			m_shiftPerpOld=- eye(num_nodes);
			v_objectiveADMM=zeros(1,obj.num_maxIterationsADMM);
			counter=0;
			
			%ADMM
			for ind_iteration = 1:obj.num_maxIterationsADMM
				
				%Y_par and Y_perp update
				m_Ypar = obj.proximalOperator( reshape(m_kroneckerDifOpPar*v_par,dim_subspace^2,dim_subspace^2)  - reshape(v_lagrangesecond,dim_subspace^2,dim_subspace^2) ...
					, obj.reg_evalMultiplicityPar/rho );
				m_Yper = obj.proximalOperator( reshape(m_kroneckerDifOpPerp*v_perp,(num_nodes-dim_subspace)^2,(num_nodes-dim_subspace)^2) ...
					- reshape(v_lagrangethird,(num_nodes-dim_subspace)^2,(num_nodes-dim_subspace)^2), obj.reg_evalMultiplicityPerp/rho );
				
				%Shift_par and shift_perp update
				v_par=m_invpar*(rho*(m_linearConspar'*v_constraint-m_linearConspar'*v_lagrangefirst+m_kroneckerDifOpPar'*v_lagrangesecond+...
					m_kroneckerDifOpPar'*m_Ypar(:)-m_linearConspar'*m_linearConsper*v_perp));
				v_perp=m_invper*(rho*(m_linearConsper'*v_constraint-m_linearConsper'*v_lagrangefirst+m_kroneckerDifOpPerp'*v_lagrangethird+...
					m_kroneckerDifOpPerp'*m_Yper(:)-m_linearConsper'*m_linearConspar*v_par));
				
				%Update of Lagrange multipliers
				v_lagrangefirst=v_lagrangefirst+(m_linearConspar*v_par+m_linearConsper*v_perp-v_constraint);
				v_lagrangesecond=v_lagrangesecond+(m_Ypar(:)-m_kroneckerDifOpPar*v_par);
				v_lagrangethird=v_lagrangethird+(m_Yper(:)-m_kroneckerDifOpPerp*v_perp);
				
				%Obtaining shift_par and shift_perp
				m_par=reshape(v_par,dim_subspace,dim_subspace);
				m_per=reshape(v_perp,num_nodes-dim_subspace,num_nodes-dim_subspace);
				m_shiftPar=U_par*m_par*U_par';
				m_shiftPerp=U_perp*m_per*U_perp';
				
				%saving the object of the optimization
				%if obj.verboseLevel >= 2
					v_objectiveADMM(1,ind_iteration) = obj.objectiveExactProjection(m_shiftPar,m_shiftPerp);
				%end
				%	 Stopping criterion
				if norm( [m_shiftPar(:);m_shiftPerp(:)] - [m_shiftParOld(:);m_shiftPerpOld(:)] )/( dim_subspace^2 + (num_nodes-dim_subspace)^2 ) < obj.toleranceADMM
					break
				end
				counter=counter+1;
				m_shiftParOld = m_shiftPar;
				m_shiftPerpOld = m_shiftPerp;
			end
			
			%output of the function
			if counter==obj.num_maxIterationsADMM % Tolerance not reached --> infeasible
				m_shift=eye(num_nodes,num_nodes);
			else
				% Project on the feasible set to ensure feasibility.
				v_x=[v_par;v_perp];
			    [m_A,v_b,m_size] = obj.feasibleSetEquations(U_par,U_perp,graph,epsilon);
                m_A=m_A(m_size+1:end,:);
                v_b=v_b(m_size+1:end,:);
				v_projected=v_x-(m_A')*((m_A*m_A')\((m_A*v_x)-v_b ));%compute projection
				m_par=reshape(v_projected(1:dim_subspace^2),dim_subspace,dim_subspace);
				m_per=reshape(v_projected(dim_subspace^2+1:end),num_nodes-dim_subspace,num_nodes-dim_subspace);
				m_shiftPar=U_par*m_par*U_par';
				m_shiftPerp=U_perp*m_per*U_perp';
				m_shift=m_shiftPar+m_shiftPerp;
                m_shift =m_shift.*(graph.m_adjacency>0);
			end
			
			% Plot objective
			if obj.verboseLevel >= 2 % extra verbose level
				
				non_zero=nnz(v_objectiveADMM);
				v_objectiveCut=v_objectiveADMM(1,1:non_zero);
				F = GFigure('m_X',1:size(v_objectiveCut,2),'m_Y', v_objectiveCut,...
					'ch_xlabel','Number of local exchanges',...
					'ch_ylabel','Objective value','c_legend',{'ADMM'});
				F.plot(3);
				pause()
				
			end
		end
		
		
		function [m_shift ,m_shiftPar,m_shiftPerp]= exactShiftSolverNoFast(obj,U_par,U_perp,graph,epsilon)
			% Quadratic solver for the scenario without nuclear norms.
			
			% input processing
			num_nodes = size(U_par,1);
			dim_subspace = size(U_par,2);
			
			
			[m_A,v_b,b_isFeasible,dim_prefeasible] = obj.feasibleSetEquations(U_par,U_perp,graph,epsilon);
			
			m_quadratic = diag([obj.reg_energyPar/(dim_subspace^2)*ones(1,dim_subspace^2),...
				obj.reg_energyPerp/(num_nodes-dim_subspace)^2*ones(1,(num_nodes-dim_subspace)^2)]);
			if obj.verboseLevel
				b_isFeasible,dim_prefeasible
				num_missingEdges = sum(sum( graph.m_adjacency == 0 ) )/2;
				fprintf('%d missing edges out of %d (%.2f %%)',num_missingEdges,(num_nodes^2 - num_nodes)/2,num_missingEdges/((num_nodes^2 - num_nodes)/2)*100);
				opt = optimoptions('quadprog','Algorithm','interior-point-convex','Display','final');
			else
				opt = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
			end
			
			v_f = quadprog(m_quadratic + 10*eps*eye(size(m_quadratic)),[],[],[],m_A,v_b,[],[],[],opt);
			
			F_par = reshape( v_f(1:dim_subspace^2) , dim_subspace , dim_subspace);
			F_perp = reshape( v_f(dim_subspace^2+1:end) , num_nodes-dim_subspace , num_nodes-dim_subspace);
			
			m_shiftPar=U_par*F_par*U_par';
			m_shiftPerp=U_perp*F_perp*U_perp';
			m_shift =  m_shiftPar+ m_shiftPerp;
		end
		
		
		function [m_shift] = approximateShiftSolver(obj,U_par,U_perp,graph,epsilon)
			%
			% m_shift is the function that minimizes the approximate
			% projection criterion described in the properties section.
			%
			% The algorithm implemented here is based on the scaled form of
			% ADMM.
			%
			% This is an initial implementation. Some performance improvements are
			% still possible.
			
			% Input processing
			num_nodes = size(U_par,1);
			dim_subspace = size(U_par,2);
			rho = obj.stepSizeADMM;
			
			% Matrix for m_shift update
			m_kroneckerDifOpPar = obj.operatorMatrixKoneckerDifference( dim_subspace );
			m_kroneckerDifOpPerp = obj.operatorMatrixKoneckerDifference( num_nodes - dim_subspace );
			m_quadraticSTerm = zeros(num_nodes^2,num_nodes^2);
			m_kroneckerDifOpFPar = m_kroneckerDifOpPar*kron(U_par',U_par');
			m_kroneckerDifOpFPerp = m_kroneckerDifOpPerp*kron(U_perp',U_perp');
			m_quadraticSTerm = m_quadraticSTerm ...
				+ rho/2* (obj.reg_evalMultiplicityPar/dim_subspace^2)^2 * (m_kroneckerDifOpFPar'*m_kroneckerDifOpFPar)...
				+ rho/2* (obj.reg_evalMultiplicityPerp/ (num_nodes-dim_subspace)^2 )^2 * (m_kroneckerDifOpFPerp'*m_kroneckerDifOpFPerp) ;
			[m_constraints,v_constraints] = obj.feasibleSetEquations(U_par,U_perp,graph,epsilon);
			m_quadraticSTerm = m_quadraticSTerm + rho/2 * (m_constraints'*m_constraints);
			U_par_par = kron(U_par,U_par);
			U_perp_perp = kron(U_perp,U_perp);
			U_par_perp = kron(U_par,U_perp);
			U_perp_par = kron(U_perp,U_par);
			m_quadraticSTerm = m_quadraticSTerm ...
				+ ( obj.reg_separation / ( (num_nodes - dim_subspace)*dim_subspace ) )...
				*.5*(U_par_perp*U_par_perp'+U_perp_par*U_perp_par')... % balancing
				+ ( obj.reg_energyPar / ( dim_subspace^2 ) )*(U_par_par*U_par_par')...
				+ ( obj.reg_energyPerp / ( num_nodes - dim_subspace)^2  )*(U_perp_perp*U_perp_perp');
			
			% Initializations
			m_shiftOld = - eye(num_nodes);
			m_Y1 = zeros( dim_subspace^2 );
			m_Y2 = zeros( (num_nodes-dim_subspace)^2 );
			m_Q11 = zeros( dim_subspace^2 );
			m_Q12 = zeros( (num_nodes-dim_subspace)^2 );
			v_q2 = zeros( length(v_constraints) , 1 );
			v_objective = NaN(1,obj.num_maxIterationsADMM);
			for ind_iteration = 1:obj.num_maxIterationsADMM
				
				if (obj.verboseLevel>5) && (ind_iteration>1)
					al_beforeSUpdate = obj.approximateProjectionAugmentedLagrangian(U_par,U_perp,graph,epsilon,m_shift,m_Y1,m_Y2,m_Q11,m_Q12,v_q2)
				end
				
				% Shift update
				v_linearSTerm = (rho/2) * (...
					(obj.reg_evalMultiplicityPar / dim_subspace^2 )* (m_kroneckerDifOpFPar'*(m_Y1(:)+m_Q11(:)) )...
					+ (obj.reg_evalMultiplicityPerp / (num_nodes-dim_subspace)^2 )*(m_kroneckerDifOpFPerp'*(m_Y2(:)+m_Q12(:)) )...
					+ m_constraints'*(v_constraints-v_q2) ...
					);
				v_shift =  m_quadraticSTerm\v_linearSTerm;
				m_shift = reshape(v_shift,num_nodes,num_nodes);
				
				if (obj.verboseLevel>5) && (ind_iteration>1)
					al_afterSUpdate = obj.approximateProjectionAugmentedLagrangian(U_par,U_perp,graph,epsilon,m_shift,m_Y1,m_Y2,m_Q11,m_Q12,v_q2)
				end
				
				% Y update
				m_FPar = U_par'*m_shift*U_par;
				m_kronDifFPar = (obj.reg_evalMultiplicityPar / dim_subspace^2 )	*(kron(m_FPar,eye(dim_subspace)) - kron(eye(dim_subspace),m_FPar  ));
				m_FPerp = U_perp'*m_shift*U_perp;
				m_kronDifFPerp = (obj.reg_evalMultiplicityPerp / (num_nodes-dim_subspace)^2 )*(kron(m_FPerp,eye(num_nodes-dim_subspace)) - kron(eye(num_nodes-dim_subspace),m_FPerp  ) );
				m_Y1 = obj.proximalOperator( m_kronDifFPar - m_Q11 , 1/rho );
				m_Y2 = obj.proximalOperator( m_kronDifFPerp- m_Q12 , 1/rho );
				
				if (obj.verboseLevel>5) && (ind_iteration>1)
					al_afterYUpdate = obj.approximateProjectionAugmentedLagrangian(U_par,U_perp,graph,epsilon,m_shift,m_Y1,m_Y2,m_Q11,m_Q12,v_q2)
				end
				
				% Update of Lagrange multipliers
				m_Q11 = m_Q11 + m_Y1 - m_kronDifFPar;
				m_Q12 = m_Q12 + m_Y2 - m_kronDifFPerp;
				v_q2 = v_q2 + m_constraints*m_shift(:) - v_constraints;
				
				if obj.verboseLevel >= 2
					v_objective(1,ind_iteration) = obj.objectiveApproximateProjection(U_par,((m_shift+m_shift')/2).*(graph.m_adjacency>0));
				end
				
				% Stopping criterion
				if norm( m_shift - m_shiftOld , 'fro' ) / num_nodes^2 < obj.toleranceADMM
					break
				end
				m_shiftOld = m_shift;
			end
			
			% Ensure the output is symmetric and topologically feasible
			m_shift = ((m_shift+m_shift')/2).*(graph.m_adjacency>0);
			% Plot objective and evals
			if obj.verboseLevel >= 2
				F(1) = GFigure('m_Y',v_objective,'ch_title','Approx. Proj. Objective -- ADMM',...
					'ch_ylabel','Objective value','ch_xlabel','Iteration index');
				m_Y = GFigure.formMatrixByPaddingRows( eig(m_FPar+m_FPar')'/2 ,  eig(m_FPerp+m_FPerp')'/2 );
				F(2) = GFigure('m_Y',m_Y,'ch_title','Eigenvalues of ADMM solution',...
					'ch_ylabel','Eigenvalues','ch_xlabel','Eigenvalue index','ch_plotType2D','stem',...
					'c_legend',{'Fpar','Fperp'},'ch_legendPosition','northwest');
				F.plot(3);
			end
			
			
		end
		
		
		function [m_shift] = approximateShiftSolverNoFast(obj,U_par,U_perp,graph)
			%
			% Quadratic solver for the scenario without nuclear norms.
			%
			
			% input processing
			num_nodes = size(U_par,1);
			dim_subspace = size(U_par,2);
			epsilon = obj.epsilon;  % will be an input arg in the future
			
			% Matrices for feasible set
			[m_A,v_b] = obj.feasibleSetEquations(U_par,U_perp,graph,epsilon);
			
			% Matrices for the objective
			m_quadratic = (obj.reg_separation/dim_subspace^2)*kron(U_par*U_par',U_perp*U_perp')...
				+ (obj.reg_energyPar/dim_subspace^2)*kron(U_par*U_par',U_par*U_par')...
				+ (obj.reg_energyPerp/((num_nodes-dim_subspace)^2))*kron(U_perp*U_perp',U_perp*U_perp');
			
			m_quadratic = m_quadratic + 100*eps*eye(size(m_quadratic,1)); % small correction for numerical errors
			
			% Solve quadratic optimization problem
			if obj.verboseLevel
				[~,~,b_isFeasible,dim_prefeasible] = obj.feasibleSetEquations(U_par,U_perp,graph,epsilon);
				b_isFeasible,dim_prefeasible
				num_missingEdges = sum(sum( graph.m_adjacency == 0 ) )/2;
				fprintf('%d missing edges out of %d (%.2f %%)',num_missingEdges,(num_nodes^2 - num_nodes)/2,num_missingEdges/((num_nodes^2 - num_nodes)/2)*100);
				opt = optimoptions('quadprog','Algorithm','interior-point-convex','Display','final');
			else
				opt = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
			end
			
			v_s = quadprog(m_quadratic,[],[],[],m_A,v_b,[],[],[],opt);
			
			if obj.verboseLevel
				fprintf('Final objective value = %g\n',v_s'*m_quadratic*v_s);
			end
			
			m_shift = reshape(v_s,num_nodes,num_nodes);
		end
		
		

		
		
		function [m_shift] = approximateShiftSolverOld(obj,m_basisSubspace,graph)
			
			%obtaining the missing edges
			[m_missingEdges] = graph.getMissingEdges();
			
			%input processing for expressing the cost function in vector form
			[U_par,U_perp] = obj.findOrthonormalBases(m_basisSubspace);
			[m_parallel]=obj.nuclearKron(U_par);% vec^-1(m_parallel*vec(m_shift)) is the another expression of nuclear_kron(U_par'*m_shift*U_par)
			[m_perp]=obj.nuclearKron(U_perp);% vec^-1(m_perp*vec(m_shift)) is the another expression of nuclear_kron(U_perp'*m_shift*U_perp)
			num_nodes=size(U_par,1);
			subspaceDimension=size(U_par,2);
			
			%input processing for the linear constraints
			m_zeroEntries=symmetric_matrix_nullspace(num_nodes);%m_zeroEntries*vec(m_shift) makes some entries of m_shift zero which belongs disconnected nodes
			ide_kron = kron(eye(num_nodes),eye(num_nodes));
			m_symmetry = ide_kron(  m_missingEdges(2,:)+num_nodes*(m_missingEdges(1,:)-1), :);%m_zeroEntries*vec(m_shift) makes m_shift symmetric
			m_feasibility=[m_symmetry;m_zeroEntries;(U_par(:))'*kron(U_par',eye(num_nodes))];%vec(U_perp)'*kron(U_perp',eye(num_nodes))];% it encompasses all constraints related to the feasibility
			v_feasibility=[zeros(size(m_symmetry,1),1);zeros(size(m_zeroEntries,1),1);subspaceDimension];%(1-obj.epsilon)*(num_nodes-subspaceDimension)];%m_feasibilit*vec(m_shift)=v_feasibility
			m_constraint_shift=[obj.reg_evalMultiplicityPar*m_parallel;zeros(subspaceDimension^2*(num_nodes-subspaceDimension)^2,num_nodes^2);...
				zeros(subspaceDimension^2*(num_nodes-subspaceDimension)^2,num_nodes^2);obj.reg_evalMultiplicityPerp*m_perp]; %it is used to expresses m_nuclearKron based on m_shift
			
			%input processing for ADMM iterations
			m_shift=zeros(num_nodes,num_nodes);
			m_shift=(m_shift+m_shift')/2;
			m_shiftOld=-m_shift;
			m_nuclearKron=[obj.reg_evalMultiplicityPar*reshape(m_parallel*m_shift(:),subspaceDimension^2,subspaceDimension^2),...
				zeros(subspaceDimension^2,(num_nodes-subspaceDimension)^2);zeros((num_nodes-subspaceDimension)^2,subspaceDimension^2)...
				,obj.reg_evalMultiplicityPerp*reshape(m_perp*m_shift(:),(num_nodes-subspaceDimension)^2,(num_nodes-subspaceDimension)^2)];
			% nuclear_norm (m_nuclearKron)=nuclear_kron(U_perp'*m_shift*U_perp)+ nuclear_kron(U_par'*m_shift*U_par)
			v_first_lagrang=zeros(size(m_nuclearKron(:)));% The first Lagrange multiplier
			v_second_lagrang=zeros(size(v_feasibility));% The second Lagrange multiplier
			v_objective=NaN(1,obj.num_iterationsADMM);
			%ADMM iteration
			for ind_numADMM=1:obj.num_iterationsADMM
				%m_shift iterations
				v_first=(obj.reg_firstAugmented/2)*(m_constraint_shift)'*m_constraint_shift...
					+(obj.reg_secondAugmented/2)*(m_feasibility'*m_feasibility)...
					+obj.reg_separation*(kron(U_par',U_perp'))'*kron(U_par',U_perp')...
					+obj.reg_energyPar*(kron(U_par',U_par'))'*kron(U_par',U_par')+obj.reg_energyPerp*(kron(U_perp',U_perp'))'*kron(U_perp',U_perp');...
					%the first part of the formulation to compute m_shift
				
				[v_lagrangkron]=obj.lagrangKron(v_first_lagrang,U_par,U_perp);% computes m_constraint_shift'*m_constraint_nuclear_inv*v_first_lagrang
				v_vec1=(U_par*U_par'*m_shift-trace(U_par*U_par'*m_shift)/subspaceDimension*eye(num_nodes))*U_par*U_par';
                v_vec2=(U_perp*U_perp'*m_shift-trace(U_perp*U_perp'*m_shift)/(num_nodes-subspaceDimension)*eye(num_nodes))*U_perp*U_perp';
				v_second=obj.reg_firstAugmented/2*(obj.reg_evalMultiplicityPar^2*2*subspaceDimension*v_vec1(:)+...
					obj.reg_evalMultiplicityPerp^2*2*(num_nodes-subspaceDimension)*v_vec2(:)+v_lagrangkron)...
					+(obj.reg_secondAugmented/2)*(m_feasibility'*v_feasibility-m_feasibility'*v_second_lagrang);%the second part of the formulation to compute m_shift
				v_shift=inv(v_first)*v_second;%vec(m_shift)
				
				%m_nuclearKron iterations
				[v_kronShift,m_kronShift]=obj.krontruncate(reshape(v_shift,num_nodes,num_nodes),U_par,U_perp,obj.reg_evalMultiplicityPar,obj.reg_evalMultiplicityPerp);...
					%computes m_constraint_nuclear_inv*m_constraint_shift*vec(m_shift)
				m_second=reshape(v_first_lagrang,(num_nodes-subspaceDimension)^2+subspaceDimension^2,(num_nodes-subspaceDimension)^2+subspaceDimension^2);%vec^-1(v_first)
				[m_leftSVD,m_shiftSVD,m_rightSVD]=svd(m_kronShift+m_second);...
					%SVD of vec^-1(inv(m_constraint_nuclear)*m_constraint_shift*v_shift)+%vec^-1(v_first) to compute the proximal operator of the nuclear norm
				v_shiftSVD=diag(m_shiftSVD);
				for ind_svd=1:length(v_shiftSVD)
					v_newshift_SVD(ind_svd)=max(v_shiftSVD(ind_svd)-inv(obj.reg_firstAugmented),0);
				end
				m_newshiftSVD=diag(v_newshift_SVD);
				m_nuclearKron=m_leftSVD*m_newshiftSVD*m_rightSVD;
				v_first_lagrang=v_first_lagrang+obj.reg_firstAugmented*(m_nuclearKron(:)-v_kronShift);
				v_second_lagrang=v_second_lagrang+obj.reg_secondAugmented*(m_feasibility*v_shift-v_feasibility);
				m_shift=reshape(v_shift,num_nodes,num_nodes);
				% Objective evaluation
				if obj.verboseLevel
					v_objective(1,ind_numADMM) = obj.objectiveApproximateProjection(m_basisSubspace,m_shift); % speed is not important in verbose mode
				end
				
				% Stopping criterion
				if norm(m_shift-m_shiftOld)<obj.toleranceADMM
					break
				end
				m_shiftOld=m_shift;
			end
			m_shift=reshape(v_shift,num_nodes,num_nodes);
			
			% Plot objective
			if obj.verboseLevel >= 2 % extra verbose level
				
				F = GFigure('m_X',1:obj.num_iterationsADMM,'m_Y', v_objective,...
					'ch_xlabel','Number of local exchanges',...
					'ch_ylabel','Objective value','c_legend',{'ADMM'});
				F.plot(3);
				pause()
				
			end
		end
		
	
	end
	
	methods(Static,Access = private)
		
		function [m_output]=subgradientNuclearKronecker(m_input)
			% The suibgradient of nucleaer_kronecker operator
			dimInput=size(m_input,1);
			[m_leftEig,m_rightEig]=eig(m_input);
			m_diagonal=kron(diag(m_rightEig),ones(dimInput,1))-kron(ones(dimInput,1),diag(m_rightEig));
			
			for ind_rows=1:dimInput^2
				if m_diagonal(ind_rows,:)>=0;
					m_diagonal(ind_rows,:)=1;
				else
					m_diagonal(ind_rows,:)=-1;
				end
			end
			for ind_first=1:dimInput
				v_firstSum(1,ind_first)=sum(m_diagonal(dimInput*(ind_first-1)+1:1:dimInput*ind_first));
			end
			for ind_second=1:dimInput
				v_secondSum(1,ind_second)=sum(m_diagonal(ind_second:dimInput:dimInput*(dimInput-1)+ind_second));
			end
			v_shift=[];
			for ind_row=1:dimInput
				for ind_column=1:dimInput
					v_shift=[v_shift ( sum(m_leftEig(ind_row,:).*m_leftEig(ind_column,:).*v_firstSum) - sum(m_leftEig(ind_row,:).*m_leftEig(ind_column,:).*v_secondSum) )];...
						% obtain the subgradient of our cost function
				end
			end
			m_output=(reshape(v_shift,dimInput,dimInput));%obtain vec^(-1)(v_shift) which is subgradient of our cost function
			
		end
		function [FF_psgd,FF_per_psgd]=projection_FF(F0,U_parallel,S0,ep,m_missingEdges,U_perp)
			r=size(U_parallel,2);
			N=size(U_parallel,1);
			[C,b]=FastDecentralizedProjectionEstimator.linear_constraint_computation(...
				U_parallel,ep,m_missingEdges,U_perp);
			x0=[F0(:);S0(:)];
			x=x0-(C')*((C*C')\((C*x0)-b));%compute projection
			FF_psgd=reshape(x(1:r^2,1),r,r);%obtain FF from the projection
			FF_per_psgd=reshape(x(r^2+1:(N-r)^2+r^2,1),N-r,N-r);%obtain S_per from the projection
		end
		
		
		function [C,b]=linear_constraint_computation(U_parallel,epsilon,m_missingEdges,U_perp)
			% rewriting the constraints of the proposed method as just one
			% linear constraint
			
			num_nodes=size(U_parallel,1);
			num_missEdges=size(m_missingEdges,2);
			subspaceDimension=size(U_parallel,2);
			assert(num_nodes >= subspaceDimension);
			[G_perp]=symmetric_matrix_nullspace(num_nodes-subspaceDimension);
			[G_parallel]=symmetric_matrix_nullspace(subspaceDimension);
			U_par_kron_U_par = kron(U_parallel,U_parallel);
			U_perp_kron_U_perp = kron(U_perp,U_perp);
			W_U_par_kron_U_par = U_par_kron_U_par(  m_missingEdges(2,:)+num_nodes*(m_missingEdges(1,:)-1)   , :);
			W_U_perp_kron_U_perp = U_perp_kron_U_perp(  m_missingEdges(2,:)+num_nodes*(m_missingEdges(1,:)-1)   , :);
			m_eye1=eye(subspaceDimension);
            m_eye2=eye(num_nodes-subspaceDimension);
			C=[W_U_par_kron_U_par W_U_perp_kron_U_perp;G_parallel zeros((subspaceDimension^2-subspaceDimension)/2,(num_nodes-subspaceDimension)^2);zeros(((num_nodes-subspaceDimension)^2-(num_nodes-subspaceDimension))/2,subspaceDimension^2) ...
				G_perp;m_eye1(:)' zeros(1,(num_nodes-subspaceDimension)^2)];
			b=[zeros(num_missEdges,1);zeros((subspaceDimension^2-subspaceDimension)/2,1);zeros(((num_nodes-subspaceDimension)^2-(num_nodes-subspaceDimension))/2,1);10*subspaceDimension]
			
			if ~isnan(epsilon)
				C = [C;zeros(1,subspaceDimension^2) m_eye2(:)'];
				b=[b;(num_nodes-subspaceDimension)*(1+epsilon)];
			end
			
			[m_topologyNullSpace] = FastDecentralizedProjectionEstimator.symmetricShiftNullSpace(num_nodes,m_missingEdges);
			m_preFeasibleA = [m_topologyNullSpace;symmetric_matrix_nullspace(num_nodes)];
			
			% form a matrix and vector to avoid ambiguity in the solution and avoid
			% trivial solution = identity
			m_trace_par = U_parallel*U_parallel';
			m_trace_perp = U_perp*U_perp';
			m_ambiguity = [ m_trace_par(:)'];
			v_ambiguity = [subspaceDimension];
			v_ambiguity = [v_ambiguity ;(1+epsilon)*(num_nodes-subspaceDimension)];
			m_ambiguity = [m_ambiguity ; m_trace_perp(:)'];
			% output computation
			m_A = [m_preFeasibleA;m_ambiguity];
			v_b = [zeros(size(m_preFeasibleA,1),1);v_ambiguity];
			
			
		end
		
		
		function [m_A] = symmetricShiftNullSpace(num_nodes,m_missingEdges)
			%
			% INPUT:
			%   num_nodes: the number of nodes
			%   m_missingEdges= contains the missing edges
			%
			% OUTPUT:
			%   m_A : M x N^2 matrix such that a symmetric matrix m_S is a
			%       valid shift matrix iff m_A * m_S(:) = 0
			%
			
			num_missingEdges = size(m_missingEdges,2);
			m_A = zeros(num_missingEdges,num_nodes^2);
			for ind_missingEdges = 1:num_missingEdges
				node_1 = m_missingEdges(1,ind_missingEdges);
				node_2 = m_missingEdges(2,ind_missingEdges);
				m_A(ind_missingEdges, num_nodes*(node_1-1) + node_2 ) = 1;
			end
			
		end
		
		
		function [m_output]=nuclearKron(m_input)
			col=1;
			dim=size(m_input,2);
			A=NaN(dim^4,dim^2);
			ide_para=eye(dim);
			for j=1:dim
				for i=1:dim
                    m_kron=kron(ide_para(:,i)*ide_para(:,j)',eye(dim))-kron(eye(dim),ide_para(:,i)*ide_para(:,j)');
					A(:,col)=m_kron(:);
					col=col+1;
				end
			end
			m_output=A*kron(m_input',m_input');
		end
		
		function [m_output]=constraintNuclear(N,r)
			T_11=[eye(r^2),zeros(r^2,(N-r)^2)];
			T_21=[zeros((N-r)^2,r^2),eye((N-r)^2,(N-r)^2)];
			T_22=[eye(r^2),zeros(r^2,(N-r)^2)];
			T_31=[eye(r^2),zeros(r^2,(N-r)^2)];
			T_32=[zeros((N-r)^2,r^2),eye((N-r)^2,(N-r)^2)];
			T_41=[zeros((N-r)^2,r^2),eye((N-r)^2,(N-r)^2)];
			T_42=[zeros((N-r)^2,r^2),eye((N-r)^2,(N-r)^2)];
			m_output=[kron(T_11,T_11);kron(T_22,T_21);kron(T_32,T_31);kron(T_42,T_41)];
		end
		
		function [v_output,m_output]=krontruncate(m_input,U_parallel,U_perp,reg_evalMultiplicityPar,reg_evalMultiplicityPerp)
			num_nodes=size(U_parallel,1);
			subspaceDim=size(U_parallel,2);
			m_first1=kron(U_parallel'*m_input*U_parallel,eye(subspaceDim));
			m_first2=kron(eye(subspaceDim),U_parallel'*m_input*U_parallel);
			m_first=reg_evalMultiplicityPar*(m_first1-m_first2);
			m_second=zeros(subspaceDim^2,(num_nodes-subspaceDim)^2);
			m_third=zeros((num_nodes-subspaceDim)^2,subspaceDim^2);
			m_forth1=kron(U_perp'*m_input*U_perp,eye(num_nodes-subspaceDim));
			m_forth2=kron(eye((num_nodes-subspaceDim)),U_perp'*m_input*U_perp);
			m_forth=reg_evalMultiplicityPerp*(m_forth1-m_forth2);
			m_output=[m_first,m_second;m_third,m_forth];
			v_output=m_output(:);
		end
		
		
		function [v_output]=lagrangKron(v_input,U_par,U_perp,reg_evalMultiplicityPar,reg_evalMultiplicityPerp)
			
			%input processing
			num_nodes=size(U_par,1);
			subspaceDim=size(U_par,2);
			T_11=[eye(subspaceDim^2),zeros(subspaceDim^2,(num_nodes-subspaceDim)^2)];
			T_41=[zeros((num_nodes-subspaceDim)^2,subspaceDim^2),eye((num_nodes-subspaceDim)^2,(num_nodes-subspaceDim)^2)];
			T_42=[zeros((num_nodes-subspaceDim)^2,subspaceDim^2),eye((num_nodes-subspaceDim)^2,(num_nodes-subspaceDim)^2)];
			m_par=NaN(subspaceDim);
			m_perp=NaN((num_nodes-subspaceDim));
			m_input=reshape(v_input,subspaceDim^2+(num_nodes-subspaceDim)^2,subspaceDim^2+(num_nodes-subspaceDim)^2);
			m_first=T_11*m_input*T_11';
			m_second=T_41*m_input*T_42';
			
			%program
			for ind_row=1:subspaceDim
				for ind_column=1:subspaceDim
					m_par(ind_row,ind_column)=trace(m_first)*(m_first(ind_row,ind_column)-m_first(ind_column,ind_row));
				end
			end
			for ind_row=1:(num_nodes-subspaceDim)
				for ind_column=1:(num_nodes-subspaceDim)
					m_perp(ind_row,ind_column)=trace(m_second)*(m_second(ind_row,ind_column)-m_second(ind_column,ind_row));
				end
            end
            m_par=U_par*m_par*U_par';
            m_per=U_perp*m_perp*U_perp';
			v_outpar=m_par(:);
			v_outperp=m_per(:);
			v_output=reg_evalMultiplicityPar*v_outpar+reg_evalMultiplicityPerp*v_outperp;
		end
		
		%% Functions for approximate
		
		function m_A = operatorMatrixKoneckerDifference( dim_input )
			% m_A is a matrix such that
			%
			% vecinv( m_A*vec(m_X)) = kron(m_X,eye(dim_input) - kron(eye(dim_input),m_X)
			%
			% for an arbitrary matrix m_X.
			%
			m_A = NaN(dim_input^4,dim_input^2);
			for ind1 = 1:dim_input
				for ind2 = 1:dim_input
					m_column = zeros( dim_input^2 , dim_input^2 );
					m_column( dim_input*(ind1-1)+1:dim_input*ind1 , dim_input*(ind2-1)+1:dim_input*ind2 ) = ...
						eye(dim_input);
					for ind_block = 1:dim_input
						m_column(dim_input*(ind_block-1)+ind2, dim_input*(ind_block-1)+ind1 ) = ...
							m_column(dim_input*(ind_block-1)+ind2, dim_input*(ind_block-1)+ind1 ) -1;
					end
					m_A(:,dim_input*(ind1-1) + ind2) = m_column(:);
				end
			end
			
		end
		
		
		function m_out = proximalOperator( m_in , thresh )
			% proximal operator of the nuclear norm evaluated at m_in and
			% with threshold thresh. See ?[cai2010singular].
			[U,S,V] = svd(m_in);
			v_svals = max(0,diag(S)-thresh);
			m_out = U*diag(v_svals)*V';
			
		end
		
	end
	
	
	methods(Static) % General utilities
		
		function [relEvalSeparation,v_evals1,v_evals2] = relativeEigenvalueSeparation(m_basisSubspace,m_shift)
			%
			% relEvalSeparation = min_distance( v_evals1 , v_evals2 )/max_distance( [v_evals1,v_evals2] )
			%
			% where v_evals1 is the vec of evals of F_par and v_evals2 is
			% the vec of evals of F_perp
			
			% Orthonormal bases
			[U_par,U_perp] = DecentralizedProjectionEstimator.findOrthonormalBases(m_basisSubspace);
			
			% Obtain eval vecs
			F_par = U_par'*m_shift*U_par;
			F_perp = U_perp'*m_shift*U_perp;
			
			v_evals1 = eig( (F_par+F_par')/2 );
			v_evals2 = eig( (F_perp+F_perp')/2 );
			
			m_dist_1_to_2 = abs(repmat(v_evals1,1,length(v_evals2)) - repmat(v_evals2',length(v_evals1),1));
			
			v_evals_1_and_2 = [v_evals1;v_evals2];
			m_dist_1_and_2 = abs(repmat(v_evals_1_and_2,1,length(v_evals_1_and_2)) - ...
				repmat(v_evals_1_and_2',length(v_evals_1_and_2),1) ); % this can be improved for performance
			
			relEvalSeparation = min(m_dist_1_to_2(:)) / max( m_dist_1_and_2(:) ) ;
			
		end
		
		
		function out = nuclearKronecker( m_in )
			
			sz = size(m_in,1);
            out = norm( svd(kron(m_in,eye(sz)) - kron(eye(sz),m_in)),1 );
			%out = norm_nuc( kron(m_in,eye(sz)) - kron(eye(sz),m_in) );
		end
		
	end
	
	
	methods % Test methods
		
		
		function val = objectiveExactProjection(obj,m_shiftPar,m_shiftPerp)
			
			dim_subspace = size(m_shiftPar,1);
			num_nodes = dim_subspace + size(m_shiftPerp,1);
			val = obj.reg_evalMultiplicityPar * obj.nuclearKronecker( m_shiftPar )  / dim_subspace^2 ...
				+ obj.reg_evalMultiplicityPerp * obj.nuclearKronecker( m_shiftPerp ) / (num_nodes - dim_subspace)^2 ...
				+ obj.reg_energyPar * norm( m_shiftPar , 'fro' )^2  / dim_subspace^2 ...
				+ obj.reg_energyPerp * norm( m_shiftPerp , 'fro' )^2 / (num_nodes - dim_subspace)^2;
			
		end
		
		
		function val = objectiveApproximateProjection(obj,m_basisSubspace,m_shift)
			
			[U_par,U_perp] = obj.findOrthonormalBases(m_basisSubspace);
			dim_subspace = size(U_par,2);
			num_nodes = size(m_shift,1);
			val = obj.objectiveExactProjection(U_par'*m_shift*U_par,U_perp'*m_shift*U_perp)+ ...
				+ obj.reg_separation * norm( U_par' * m_shift * U_perp ,'fro')^2  / ( (num_nodes - dim_subspace) * dim_subspace );
			
		end
		
		
		function val = approximateProjectionAugmentedLagrangian(obj,U_par,U_perp,graph,epsilon,m_shift,m_Y1,m_Y2,m_Q11,m_Q12,v_q2)
			% Augmented Lagrangian method to test ADMM for approximate projections
			
			num_nodes = size(U_par,1);
			dim_subspace = size(U_par,2);
			rho = obj.stepSizeADMM;
			[m_constraints,v_constraints] = obj.feasibleSetEquations(U_par,U_perp,graph,epsilon);
			
			%  Eval multiplicity terms
			val = norm_nuc(m_Y1) + norm_nuc(m_Y2);
			
			% Separation term
			val = val + (obj.reg_separation  / ( (num_nodes - dim_subspace)* dim_subspace ) )...
				* 0.5*(norm( U_par' * m_shift * U_perp ,'fro')^2  + norm( U_perp' * m_shift * U_par ,'fro')^2  );
			
			% Energy terms
			val = val + (obj.reg_energyPar/ dim_subspace^2) * norm( U_par' * m_shift * U_par ,'fro' )^2 ...
				+ (obj.reg_energyPerp / (num_nodes - dim_subspace)^2) * norm( U_perp' * m_shift * U_perp ,'fro' )^2;
			
			% Augmented term for Y1 and Y2
			m_FPar = U_par'*m_shift*U_par;
			m_FPerp = U_perp'*m_shift*U_perp;
			val = val + (rho/2)*norm( m_Y1 - (obj.reg_evalMultiplicityPar / dim_subspace^2 ) ...
				* (kron(m_FPar,eye(dim_subspace))-kron(eye(dim_subspace),m_FPar) ) + m_Q11   ,'fro')^2 ...
				+ (rho/2)*norm( m_Y2 - (obj.reg_evalMultiplicityPerp / (num_nodes-dim_subspace)^2 ) ...
				* (kron(m_FPerp,eye(num_nodes-dim_subspace))-kron(eye(num_nodes -dim_subspace),m_FPerp) ) + m_Q12   ,'fro')^2;
			
			% Augmented term for v_q2
			val = val + (rho/2)*norm( m_constraints*m_shift(:) - v_constraints + v_q2)^2;
			
		end
		
		
		%% Test methods under construction
		
		function F = plotSeparationVsWeightParallel(obj,m_basisSubspace,graph,v_weightParallel,v_epsilon)
			%
			% Retuns a GFigure for analyzing the distance between evals.
			% The behavior depends on whether v_epsilon is a vector or a
			% scalar.
			
			% input processing
			num_nodes = size(m_basisSubspace,1);
			dim_subspace = size(m_basisSubspace,2);
			
			for ind_epsilon = length(v_epsilon):-1:1
				for ind_weight = 1:length(v_weightParallel)
					[F_par,F_perp]= obj.exactShiftEvalRegularizedGivenWeight(m_basisSubspace,graph,v_weightParallel(ind_weight),v_epsilon(ind_epsilon));
					Eigs_par(:,ind_weight) = sort(eig((F_par+F_par')/2));
					Eigs_perp(:,ind_weight) = sort(eig((F_perp+F_perp')/2));
					relEvalSeparation(ind_epsilon,ind_weight) = obj.relativeEigenvalueSeparation(Eigs_par(:,ind_weight),Eigs_perp(:,ind_weight));
				end
				
				%debug
				vpar=(Eigs_par(:,ind_weight));
				vper=(Eigs_perp(:,ind_weight));
				mx = max([vpar',vper']);
				mi = min([vpar',vper']);
				%vpar = (vpar - (mx+mi)/2)/( mx  - (mx+mi)/2 );
				%vper = (vper - (mx+mi)/2)/( mx  - (mx+mi)/2 );
				F_perp_norm = sort(eig(F_perp - (mx+mi)/2*eye(size(F_perp)))/( mx  - (mx+mi)/2 ))
				%[vper F_perp]
				%vper - sort(eig((F_perp+F_perp')/2))
				
				
				%vper
				
				c_legend{ind_epsilon} = sprintf('\\epsilon = %g',v_epsilon(ind_epsilon));
			end
			
			m_multiplot(1,1) = GFigure('m_X',v_weightParallel,'m_Y',relEvalSeparation,...
				'ch_title','Relative Eigenvalue Separation','ch_xlabel','weightParallel','c_legend',c_legend);
			if length(v_epsilon) == 1
				m_color = [zeros(dim_subspace,3);ones(num_nodes-dim_subspace,3)*diag([1,0,0])];
				m_multiplot(2,1) = GFigure('m_X',v_weightParallel,'m_Y',[Eigs_par;Eigs_perp],'ch_title','Black: evals(Fpar); Red: evals(Fperp)','m_colorset',m_color,...
					'c_styles',[repmat({'-'},1,dim_subspace),repmat({'--'},1,num_nodes-dim_subspace)],...
					'colorPeriod',num_nodes,'ch_xlabel','weightParallel');
			end
			F = GFigure('m_multiplot',m_multiplot);
			
			if 1
				
			end
			
			
			
		end
		
		
		function [m_S,m_projection] = exactShiftTuneEpsilon(obj,m_basisSubspace,graph)
			
			% input processing
			num_nodes = size(m_basisSubspace,1);
			
			% orthonormal bases
			[U_parallel,U_perp,m_projection] = obj.findOrthonormalBases(m_basisSubspace);
			
			% Tune epsilon
			v_epsilon = obj.epsilon*2.^(1:obj.maxParameterValues);
			for ind_epsilon = 1:length(v_epsilon)
				
				% obtain optimum shift for epsilon
				[m_S]= obj.exactShiftSolver(U_parallel,U_perp,graph,v_epsilon(ind_epsilon));
				
				% error full order filter
				[~,normalized_error] = obj.findCoefficientsForOrder(m_S,...
					m_projection,num_nodes-1);
				if obj.verboseLevel
					fprintf('epsilon = %g, normalized_error = %g\n',v_epsilon(ind_epsilon),normalized_error)
				end
				
				% stopping condition
				if (normalized_error/(num_nodes^2) < obj.max_errorBetweenProjectorAndFilter)
					if obj.verboseLevel
						fprintf('normalized error = %f\n',normalized_error)
					end
					return
				end
				
			end
			
			% Target error has not be achieved --> we can still return the best
			% shift so far.
			obj.printShiftMatrixAnalysis(m_S,m_basisSubspace,graph,v_epsilon(1));
			F = obj.plotEpsilonPath(m_basisSubspace,graph,v_epsilon);
			F.plot(3);
			disp('no suitable epsilon found')
			keyboard
			
		end
		
		
		function F = printShiftMatrixAnalysis(obj,m_S,m_basisSubspace,graph,epsilon)
			% This function displays useful information about the shift
			% matrix and returns a GFigure object with the evals.
			
			% feasibility info
			[~,~,b_isFeasible,dim_prefeasible] = feasibleSetEquations(obj,m_basisSubspace,graph,epsilon);
			fprintf('For epsilon = %g, Opt. problem is feasible = %d, Dim. prefeasible shift subspace = %d\n',epsilon,b_isFeasible,dim_prefeasible);
			
			% eigenvalue info
			[U_par,U_perp,m_projection] = obj.findOrthonormalBases(m_basisSubspace);
			
			F_par = U_par' * m_S * U_par;
			F_perp = U_perp' * m_S * U_perp;
			
			Eigs_par = sort(eig((F_par+F_par')/2))';
			Eigs_perp = sort(eig((F_perp+F_perp')/2))';
			
			F = GFigure('m_X',GFigure.formMatrixByPaddingRows(Eigs_par,Eigs_perp),...
				'm_Y',GFigure.formMatrixByPaddingRows(ones(1,length(Eigs_par)),.5*ones(1,length(Eigs_perp))),...
				'c_legend',{'Eigs Fpar','Eigs Fperp'},'ch_plotType2D','stem');
			
		end
		
		
		function F = plotEpsilonPath(obj,m_basisSubspace,graph,v_epsilon)
			%
			% Returns a GFigure with a plot of the evals of F_par and F_per
			% for the values of epsilon in v_epsilon
			
			[U_par,U_perp,m_projection] = obj.findOrthonormalBases(m_basisSubspace);
			
			dim_subspace = size(U_par,2);
			num_nodes = size(U_par,1);
			Eigs_par = NaN(dim_subspace,length(v_epsilon));
			Eigs_perp = NaN(num_nodes - dim_subspace,length(v_epsilon));
			
			for ind_epsilon = 1:length(v_epsilon)
				
				switch obj.str_criterion
					case('exact')
						[m_S]=	obj.exactShiftSolver(U_par,U_perp,graph,v_epsilon(ind_epsilon));
					case('approximate')
						[m_S]=	obj.approximateShiftSolver(U_par,U_perp,graph,v_epsilon(ind_epsilon));
					case('auto')
						error('not implemented')
				end
				
				F_par = U_par' * m_S * U_par;
				F_perp = U_perp' * m_S * U_perp;
				
				Eigs_par(:,ind_epsilon) = sort(eig((F_par+F_par')/2));
				Eigs_perp(:,ind_epsilon) = sort(eig((F_perp+F_perp')/2));
			end
			m_color = [zeros(dim_subspace,3);ones(num_nodes-dim_subspace,3)*diag([1,0,0])];
			F = GFigure('m_X',v_epsilon,'m_Y',[Eigs_par;Eigs_perp],'ch_title','Black: evals(Fpar); Red: evals(Fperp)','m_colorset',m_color,...
				'c_styles',[repmat({'-'},1,dim_subspace),repmat({'--'},1,num_nodes-dim_subspace)],...
				'colorPeriod',num_nodes);
			
			if nargout <1
				F.plot(2);
			end
			
			
			
		end
		
		
	end
	
end