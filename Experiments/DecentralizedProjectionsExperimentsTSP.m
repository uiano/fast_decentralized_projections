%
%  Template for experiment files. Provide here general information on your
%  simulation study, e.g.: This file contains experiments for
%  the TSP paper.
%
%

classdef DecentralizedProjectionsExperimentsTSP < ExperimentFunctionSet
	
	
	properties
		% You may define here constants that you will use in the
		% experiments of this file
		
	end
	
	
	methods
		
		%% ----------------------------------------------------------------
		%  SIMULATIONS for the paper of Fast Decentralized method
		%% ----------------------------------------------------------------
		
		%% Figure 1
		% Monte Carlo simulation which represents the error between
		% designed filter and the projection matrix vs number of
		% local exchanges.
		% - Graph: Erdos-Renyi random graph
		% - Comparison of three estimators (fast, gossiping, Rank one)
		function F = experiment_3001(obj,niter)
			% Set estimation problem parameters
			num_nodes =30;
			v_missingEdgeProb = [.7,.8];
			dim_subspace =1;
			num_localExchanges = num_nodes;
			counter=0;
			% Graph generator
			graphGen = ErdosRenyiGraphGenerator();
			graphGen.num_nodes = num_nodes ;
			graphGen.prob_selfLoop = 1;
			graphGen.b_symmetric=1;
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp = .9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon=.1;
			% Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			GossipingEstimator.b_verbose=1;
			% Construct an estimator for rank-one method
			RankOneEstimator = RankOneDecentralizedProjectionEstimator();
			RankOneEstimator.b_diag=1;
			% running
			for ind=1:length(v_missingEdgeProb)
				graphGen.prob_edge = 1-v_missingEdgeProb(ind);
				% Run the simulation
				[ v_err_proposed] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges , FastEstimator);
				[ v_err_gossiping] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges,GossipingEstimator);
				[ v_err_rankone] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges,RankOneEstimator);
				m_err_proposed(ind,:)=v_err_proposed;
				m_err_gossiping(ind,:)=v_err_gossiping;
				m_err_rankone(ind,:)=v_err_rankone;
				m_Y(ind+counter,:)= m_err_gossiping(ind,:);
				m_Y(ind+counter+1,:)= m_err_proposed(ind,:);
				m_Y(ind+counter+2,:)= m_err_rankone(ind,:);
				counter=counter+2;
			end
			% Plot the results
			F = GFigure('m_X',1:num_localExchanges,'m_Y', m_Y,'ch_xlabel','Number of local exchanges',...
				'ch_ylabel','NMPE(H_l)','c_legend',{'Gossiping, p_{miss}=.7','Exact Projection, p_{miss}=.7', 'Rank-one Projection,p_{miss}=.7', ...
				'Gossiping, p_{miss}=.8','Exact Projection, p_{miss}=.8', 'Rank-one Projection,p_{miss}=.8'});
		end
		%% Figure 2
		% Monte Carlo simulation which represents NMPE vs number of
		% local exchanges.
		% - Graph: Erdos-Renyi random graph
		% - Comparison of two estimators (fast, gossiping)
		function F = experiment_3002(obj,niter)
			
			% Set estimation problem parameters
			num_nodes =20;
			missingEdgeProb = .6;
			v_dim_subspace =[2,3,5];
			num_localExchanges = num_nodes;
			m_err_proposed=zeros(size(v_dim_subspace,2),num_nodes);
			m_err_gossiping=zeros(size(v_dim_subspace,2),num_nodes);
			counter=0;
			% Graph generator
			graphGen = ErdosRenyiGraphGenerator();
			graphGen.prob_edge = 1-missingEdgeProb;
			graphGen.num_nodes = num_nodes ;
			graphGen.prob_selfLoop = 1;
			graphGen.b_symmetric=1;
			
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp = .9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon=.1;
			
			% Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			GossipingEstimator.b_verbose=1;
			for ind=1:size(v_dim_subspace ,2)
				% Run the simulation
				[ m_err_proposed(ind,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					v_dim_subspace(ind) , num_localExchanges , FastEstimator);
				[ m_err_gossiping(ind,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					v_dim_subspace(ind) , num_localExchanges,GossipingEstimator);
				m_Y(ind+counter,:)= m_err_gossiping(ind,:);
				m_Y(ind+counter+1,:)= m_err_proposed(ind,:);
				counter=counter+1;
			end
			% Plot the results
			F = GFigure('m_X',1:num_localExchanges,'m_Y',m_Y,'ch_xlabel','Number of local exchanges',...
				'ch_ylabel','NMPE(H_l)','c_legend',{'Gossiping, r=2','Exact Projection', 'Gossiping, r=3','Exact Projection' ,'Gossiping, r=5','Exact Projection'});
		end
		%% Figure 3
		% Monte Carlo simulation which represents the error between
		% designed filter and the projection matrix vs number of
		% local exchanges.
		% - Graph: WSN random graph
		% - Comparison of two estimators (fast, gossiping)
		function F = experiment_3003(obj,niter)
			% Set estimation problem parameters
			v_num_nodes =[15,25,30];
			v_dim_subspace =[2,3,4];
			counter=0;
			threshold=.6;
			% Graph generator
			graphGen =WSNGraphGenerator();
			graphGen.threshold=threshold;
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp = .9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon=.1;
			% Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			GossipingEstimator.b_verbose=1;
			% running
			for ind=1:length(v_num_nodes)
				graphGen.num_nodes = v_num_nodes(ind);
				% Run the simulation
				[ m_err_proposed(ind,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					v_dim_subspace(ind) , v_num_nodes(end) , FastEstimator);
				[ m_err_gossiping(ind,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					v_dim_subspace(ind) , v_num_nodes(end),GossipingEstimator);
				m_Y(ind+counter,:)= m_err_gossiping(ind,:);
				m_Y(ind+counter+1,:)= m_err_proposed(ind,:);
				counter=counter+1;
			end
			
			% Plot the results
			F = GFigure('m_X',1:v_num_nodes(end),'m_Y', m_Y,'ch_xlabel','Number of local exchanges',...
				'ch_ylabel','NMPE(H_l)','c_legend',{'Gossiping, N=15, r=2','Exact Projection, N=15, r=2', ...
				'Gossiping, N=25, r=3','Exact Projection, , N=25, r=3' ,'Gossiping, N=30, r=4','Exact Projection,N=30, r=4'});
			
		end
		%% Figure 4
		% Monte Carlo simulation which represents NMSE vs number of
		% local exchanges.
		% - Graph: Erdos-Renyi random graph
		% - Comparison of four estimators (fast,gossiping,DGD,and DLMS)
		function F = experiment_3004(obj,niter)
			
			% Set estimation problem parameters
			num_nodes =20;
			missingEdgeProb = .6;
			dim_subspace =3;
			snr = 5;
			num_localExchanges = 2000;
			
			% Graph generator
			graphGen = ErdosRenyiGraphGenerator();
			graphGen.prob_edge = 1-missingEdgeProb;
			graphGen.num_nodes = num_nodes ;
			graphGen.prob_selfLoop = 1;
			graphGen.b_symmetric=1;
			
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp =.9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon=.1;
			% Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			
			%Construct an estimator for DLMS method
			DLMSEstimator=DLMSEDecentralizedProjectionEstimator();
			DLMSEstimator.stepSizeDLMS=5e-2;
			DLMSEstimator.lagrangeDLMS=1e0;
			
			%Construct an estimator for DGD method
			DGDEstimator=DGDDecentralizedProjectionEstimator();
			DGDEstimator.stepSizeDGD=1e-2;
			
			% Run the simulation
			[ NMSE_Fast]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , dim_subspace, snr , num_nodes , FastEstimator);
			[ NMSE_gossip]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , dim_subspace, snr ,  num_nodes , GossipingEstimator);
			[ NMSE_DLMS]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , dim_subspace, snr , num_localExchanges , DLMSEstimator);
			[ NMSE_DGD]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , dim_subspace, snr , num_localExchanges , DGDEstimator);
			NMSE_Fast=[NMSE_Fast NMSE_Fast(end)*ones(1,num_localExchanges-num_nodes)];
			% Plot the results
			plot(1:num_nodes,NMSE_Fast(1:num_nodes),'r','LineWidth',2)
			grid on
			hold on
			plot(1:num_nodes,NMSE_gossip,'k','LineWidth',2)
			plot(1:num_nodes,NMSE_DLMS(1:num_nodes),'b','LineWidth',2)
			plot(1:num_nodes,NMSE_DGD(1:num_nodes),'g','LineWidth',2)
			legend('Exact Projection','Gossiping','DLMS' ,'DGD');
			xlabel('Number of local exchanges','FontSize',12);
			ylabel('NMSE','FontSize',12);
			axes('position',[0.38 0.42 0.28 0.30])
			box on
			hold on
			plot(1:num_localExchanges,NMSE_DLMS,'b');
			plot(1:num_localExchanges,NMSE_DGD,'g');
			plot(1:num_localExchanges,NMSE_Fast,'r');
			ax = gca;
			ax.YAxis.Exponent = -1;
			axis tight
			box off
			F=[]
		end
		%% Figure 5
		% Monte Carlo simulation which represents NMSE vs number of
		% local exchanges.
		% - Graph: Erdos-Renyi random graph
		% - Comparison of four estimators (fast,gossiping,DGD,and DLMS)
		function F = experiment_3005(obj,niter)
			
			% Set estimation problem parameters
			num_nodes =30;
			missingEdgeProb = .7;
			dim_subspace =5;
			snr = 5;
			num_localExchanges = 2000;
			% Graph generator
			graphGen = ErdosRenyiGraphGenerator();
			graphGen.prob_edge = 1-missingEdgeProb;
			graphGen.num_nodes = num_nodes ;
			graphGen.prob_selfLoop = 1;
			graphGen.b_symmetric=1;
			
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp =.9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon=.1;
			% Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			
			%Construct an estimator for DLMS method
			DLMSEstimator=DLMSEDecentralizedProjectionEstimator();
			DLMSEstimator.stepSizeDLMS=5e-2;
			DLMSEstimator.lagrangeDLMS=1e0;
			
			%Construct an estimator for DGD method
			DGDEstimator=DGDDecentralizedProjectionEstimator();
			DGDEstimator.stepSizeDGD=1e-2;
			
			% Run the simulation
			[ NMSE_Fast]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , dim_subspace, snr , num_nodes , FastEstimator);
			[ NMSE_gossip]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , dim_subspace, snr , num_nodes , GossipingEstimator);
			[ NMSE_DLMS]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , dim_subspace, snr , num_localExchanges , DLMSEstimator);
			[ NMSE_DGD]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , dim_subspace, snr , num_localExchanges , DGDEstimator);
			NMSE_Fast=[NMSE_Fast NMSE_Fast(end)*ones(1,num_localExchanges-num_nodes)];
			% Plot the results
			plot(1:num_nodes,NMSE_Fast(1:num_nodes),'r','LineWidth',2)
			grid on
			hold on
			plot(1:num_nodes,NMSE_gossip,'k','LineWidth',2)
			plot(1:num_nodes,NMSE_DLMS(1:num_nodes),'b','LineWidth',2)
			plot(1:num_nodes,NMSE_DGD(1:num_nodes),'g','LineWidth',2)
			legend('Exact Projection','Gossiping','DLMS' ,'DGD');
			xlabel('Number of local exchanges','FontSize',12);
			ylabel('NMSE','FontSize',12);
			axes('position',[0.36 0.40 0.28 0.30])
			box on
			hold on
			plot(1:num_localExchanges,NMSE_DLMS,'b');
			plot(1:num_localExchanges,NMSE_DGD,'g');
			plot(1:num_localExchanges,NMSE_Fast,'r');
			ax = gca;
			ax.YAxis.Exponent = -1;
			axis tight
			box off
			F=[]
		end
		%% Figure 6
		% Monte Carlo simulation which represents NMSE vs the subspace
		% dimension.
		% - Graph: Erdos-Renyi random graph
		% - Comparison of four estimators (fast,gossiping,DGD,and DLMS)
		function F = experiment_3006(obj,niter)
			
			% Set estimation problem parameters
			v_num_nodes =[25,40];
			v_missingEdgeProb =.7;
			v_dim_subspace =[4,5,6,7];
			snr = 5;
			num_localExchanges =20;
			counter=0;
			% Graph generator
			graphGen = ErdosRenyiGraphGenerator();
			graphGen.prob_selfLoop = 1;
			graphGen.b_symmetric=1;
			graphGen.prob_edge = 1-v_missingEdgeProb;
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp =.9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon=.1;
			
			% Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			
			%Construct an estimator for DLMS method
			DLMSEstimator=DLMSEDecentralizedProjectionEstimator();
			DLMSEstimator.stepSizeDLMS=1e0;
			DLMSEstimator.lagrangeDLMS=1e-2;
			
			%Construct an estimator for DGD method
			DGDEstimator=DGDDecentralizedProjectionEstimator();
			DGDEstimator.stepSizeDGD=.1;
			for ind_num=1:length(v_num_nodes)
				% Graph generator
				graphGen.num_nodes = v_num_nodes(ind_num) ;
				parfor ind=1:size(v_dim_subspace,2)
					% Run the simulation
					[ NMSE_Fast]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , v_dim_subspace(ind), snr , num_localExchanges , FastEstimator);
					[ NMSE_gossip]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , v_dim_subspace(ind), snr , num_localExchanges , GossipingEstimator);
					[ NMSE_DLMS]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , v_dim_subspace(ind), snr , num_localExchanges , DLMSEstimator);
					[ NMSE_DGD]= EstimationNMSEvsNumLocalExchanges( niter , graphGen , v_dim_subspace(ind), snr , num_localExchanges , DGDEstimator);
					v_NMSE_Fast(ind_num,ind)=NMSE_Fast(end);
					v_NMSE_gossip(ind_num,ind)=NMSE_gossip(end);
					v_NMSE_DLMS(ind_num,ind)=NMSE_DLMS(end);
					v_NMSE_DGD(ind_num,ind)=NMSE_DGD(end);
				end
			end
			for ind_col=1:size(v_NMSE_Fast,1)
				m_Y(ind_col+counter,:)= v_NMSE_gossip(ind_col,:);
				m_Y(ind_col+counter+1,:)= v_NMSE_Fast(ind_col,:);
				m_Y(ind_col+counter+2,:)= v_NMSE_DLMS(ind_col,:);
				m_Y(ind_col+counter+3,:)= v_NMSE_DGD(ind_col,:);
				counter=counter+3;
			end
			% Plot the results
			F = GFigure('m_X',v_dim_subspace,'m_Y', m_Y,'ch_xlabel','The subspace dimension',...
				'ch_ylabel','NMSE','c_legend',{'Gossiping, N=25','Exact Projection, N=25',...
				'DLMS, N=25','DGD, N=25','Gossiping, N=40','Exact Projection, N=40','DLMS, N=40','DGD, N=40'});
		end
		%% Figure 7
		% Simulation to compute the number of different eigenvalues
		%  and plot the number of distinct eigenvalues versus the
		%  subspace dimension
		function F = experiment_3007(obj,niter)
			
			% Set estimation problem parameter
			v_num_nodes = [20,30,40];
			v_missingEdgeProb = [.4,.5,.6,.7,.8];
			v_dim_subspace =[2,3,4] ;
			evalTol=5e-3;
			v_num_evals=NaN(1,length(v_missingEdgeProb));
			for ind_num=1:size(v_num_nodes,2)
				for ind_miss=1:length(v_missingEdgeProb)
					% Graph generator
					graphGen = ErdosRenyiGraphGenerator();
					graphGen.prob_edge = 1-v_missingEdgeProb(ind_miss);
					graphGen.num_nodes = v_num_nodes(ind_num) ;
					graphGen.prob_selfLoop = 1;
					graphGen.b_symmetric=1;
					
					% Construct an estimator for fast decentralized method and set its parameters
					FastEstimator = FastDecentralizedProjectionEstimator();
					FastEstimator.reg_energyPar = 0;
					FastEstimator.reg_energyPerp = 1;
					FastEstimator.reg_evalMultiplicityPar = 0;
					FastEstimator.reg_evalMultiplicityPerp =1;
					FastEstimator.reg_evalMultiplicityPar = .1;
					FastEstimator.reg_evalMultiplicityPerp =.9;
					FastEstimator.num_maxIterationsADMM = 1000;
					FastEstimator.verboseLevel = 1;
					FastEstimator.stepSizeADMM = .1;
					FastEstimator.epsilon=.1;
					
					[v_num_evals(ind_num,ind_miss)] = MeanNumDistinctEvals( niter , graphGen , v_dim_subspace(ind_num), evalTol, FastEstimator);
				end
				m_Y(ind_num,:)=v_num_evals(ind_num,:);
			end
			% Plot the results
			F = GFigure('m_X',v_missingEdgeProb,'m_Y', m_Y,'ch_xlabel','p_{miss}',...
				'ch_ylabel','MNDE','c_legend',{'Exact Projection, N=20, r=2', 'Exact Projection, N=30, r=3', 'Exact Projection, N=40, r=4'});
		end
		%% Figure 8
		% Simulation to compute how many iterations we
		% need to to meet a specific value of NMPE
		function F = experiment_3008(obj,niter)
			
			% Set estimation problem parameter
			v_num_nodes = [25,30,40];
			v_dim_subspace =[3,4,5];
			v_threshold=.4:.05:.6;
			error_thresh=0.1;
			counter=0;
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = 0;
			FastEstimator.reg_evalMultiplicityPerp =1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp =.9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon=.1;
			
			% Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			GossipingEstimator.b_verbose=1;
			for ind_num=1:length(v_num_nodes)
				for ind=1:length(v_threshold)
					for ind_niter=1:niter
						% Graph generator
						graphGen =WSNGraphGenerator();
						graphGen.num_nodes = v_num_nodes(ind_num) ;
						graphGen.threshold=v_threshold(ind);
						num_localExchanges=v_num_nodes(1,ind_num);
						% Finding error
						[ v_err_proposed] = FilterNMSEvsNumLocalExchanges( 1, graphGen ,...
							v_dim_subspace(ind_num) , num_localExchanges , FastEstimator);
						[ v_err_gossiping] = FilterNMSEvsNumLocalExchanges( 1, graphGen ,...
							v_dim_subspace(ind_num) , num_localExchanges,GossipingEstimator);
						[b_proposed(ind_niter)]=find_entry(v_err_proposed,error_thresh)
						[b_gossiping(ind_niter)]=find_entry(v_err_gossiping,error_thresh)
					end
					v_containPro(ind_num,ind)=mean(b_proposed);
					v_containGos(ind_num,ind)=mean(b_gossiping);
					b_proposed=[];
					b_gossiping=[];
				end
				m_Y(ind_num+counter,:)=v_containPro(ind_num,:);
				m_Y(ind_num+counter+1,:)=v_containGos(ind_num,:);
				counter=counter+1;
			end
			% Plot figure
			F = GFigure('m_X',v_threshold,'m_Y', m_Y,'ch_xlabel','Average Number of needed local exchanges',...
				'ch_ylabel','d_{max}','c_legend',{'Exact Projection, N=25, r=3','Gossiping, N=25, r=3','Exact Projection, N=30, r=4','Gossiping, N=30, r=4'...
				, 'Exact Projection, N=40, r=5','Gossiping, N=40, r=5'});
			function [s_err]=find_entry(v_err,s_thresh)
				v_ind=find(v_err<=s_thresh);
				if numel(v_ind)>=1
					s_err=v_ind(1,1);
				else
					s_err=50;
				end
			end
			
		end
		%% Figure 9
		% Simulation to compute NMPE versus p_miss
		%  to compare the approximate and exact projection methods
		function F = experiment_3009(obj,niter)
			
			% Set estimation problem parameter
			num_nodes =20;
			dim_subspace =3;
			v_missingEdgeProb=.5:.05:.9;
			num_localExchanges=num_nodes;
			% Graph generator
			graphGen = ErdosRenyiGraphGenerator();
			graphGen.num_nodes = num_nodes;
			graphGen.prob_selfLoop = 1;
			graphGen.b_symmetric=1;
			% Construct an estimator for fast decentralized method and set its parameter
			ApproximateEstimator = FastDecentralizedProjectionEstimator();
			ApproximateEstimator.reg_energyPar = 0;
			ApproximateEstimator.reg_energyPerp = 1;
			ApproximateEstimator.reg_evalMultiplicityPar = 0;
			ApproximateEstimator.reg_evalMultiplicityPerp =1;
			ApproximateEstimator.reg_evalMultiplicityPar = .1;
			ApproximateEstimator.reg_evalMultiplicityPerp =.9;
			ApproximateEstimator.num_maxIterationsADMM = 1000;
			ApproximateEstimator.verboseLevel = 1;
			ApproximateEstimator.stepSizeADMM = .1;
			ApproximateEstimator.epsilon=.2;
			ApproximateEstimator.reg_separation = 10;
			ApproximateEstimator.str_criterion = 'approximate';
			%Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			GossipingEstimator.b_verbose=1;
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = 0;
			FastEstimator.reg_evalMultiplicityPerp =1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp =.9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon = .2;
			% algorithm
			for ind_num=1:length(v_missingEdgeProb)
				graphGen.prob_edge = 1-v_missingEdgeProb(ind_num);
				% Computing NMPE
				[ m_err_Approx(ind_num,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace, num_localExchanges , ApproximateEstimator);
				[ m_err_Gossip(ind_num,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges , GossipingEstimator);
				[ m_err_Fast(ind_num,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges , FastEstimator);
				NMPE_Approx(ind_num)= m_err_Approx(ind_num,end);
				NMPE_Gossip(ind_num)= m_err_Gossip(ind_num,end);
				NMPE_Fast(ind_num)= m_err_Fast(ind_num,end);
			end
			% Plot the results
			F = GFigure('m_X',v_missingEdgeProb,'m_Y', [ NMPE_Gossip;NMPE_Fast;NMPE_Approx],...
				'ch_xlabel','p_{miss}','ch_ylabel','NMSE','c_legend',{'Gossiping','Exact Projection','Approximate Projection'});
		end
		%% Figure 10
		% Simulation to compute NMPE versus the
		% number of local exchanges to compare the
		% approximate and exact projection methods
		function F = experiment_3010(obj,niter)
			
			% Set estimation problem parameter
			num_nodes = 25;
			dim_subspace =3;
			v_missingEdgeProb=[.85,.9];
			num_localExchanges=num_nodes;
			counter=0;
			% Graph generator
			graphGen = ErdosRenyiGraphGenerator();
			graphGen.num_nodes = num_nodes;
			graphGen.prob_selfLoop = 1;
			graphGen.b_symmetric=1;
			% Construct an estimator for fast decentralized method and set its parameter
			ApproximateEstimator = FastDecentralizedProjectionEstimator();
			ApproximateEstimator.reg_energyPar = 0;
			ApproximateEstimator.reg_energyPerp = 1;
			ApproximateEstimator.reg_evalMultiplicityPar = 0;
			ApproximateEstimator.reg_evalMultiplicityPerp =1;
			ApproximateEstimator.reg_evalMultiplicityPar = .1;
			ApproximateEstimator.reg_evalMultiplicityPerp =.9;
			ApproximateEstimator.num_maxIterationsADMM = 1000;
			ApproximateEstimator.verboseLevel = 1;
			ApproximateEstimator.stepSizeADMM = .1;
			ApproximateEstimator.epsilon=.2;
			ApproximateEstimator.reg_separation = 10;
			ApproximateEstimator.str_criterion = 'approximate';
			%Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			GossipingEstimator.b_verbose=1;
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = 0;
			FastEstimator.reg_evalMultiplicityPerp =1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp =.9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.epsilon=.2;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			% algorithm
			for ind_num=1:length(v_missingEdgeProb)
				graphGen.prob_edge = 1-v_missingEdgeProb(ind_num);
				% Computing NMPE
				[m_Y(ind_num+counter,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges , FastEstimator);
				[ m_Y(ind_num+counter+1,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace, num_localExchanges , ApproximateEstimator);
				[ m_Y(ind_num+counter+2,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges , GossipingEstimator);
				counter=counter+2;
			end
			% Plot the results
			F = GFigure('m_X',1:num_localExchanges,'m_Y', m_Y,'ch_xlabel','Number of local exchanges',...
				'ch_ylabel','NMPE(H_l)','c_legend',{'Exact Projection,p_{miss}=.85','Approximate Projection,p_{miss}=.85','Gossiping,p_{miss}=.85',...
				'Exact Projection, p_{miss}=.9' ,'Approximate Projection, p_{miss}=.9','Gossiping, p_{miss}=.9'});
		end
		%% Figure 11
		% Simulation to compute NMPE versus the
		% number of local exchanges to compare the
		% approximate and exact projection methods
		function F = experiment_3011(obj,niter)
			
			% Set estimation problem parameter
			num_nodes = 35;
			dim_subspace =5;
			v_threshold=[.25,.35];
			counter=0;
			% Graph generator
			graphGen =WSNGraphGenerator();
			graphGen.num_nodes = num_nodes ;
			num_localExchanges=num_nodes;
			% Construct an estimator for fast decentralized method and set its parameter
			ApproximateEstimator = FastDecentralizedProjectionEstimator();
			ApproximateEstimator.reg_energyPar = 0;
			ApproximateEstimator.reg_energyPerp = 1;
			ApproximateEstimator.reg_evalMultiplicityPar = 0;
			ApproximateEstimator.reg_evalMultiplicityPerp =1;
			ApproximateEstimator.reg_evalMultiplicityPar = .1;
			ApproximateEstimator.reg_evalMultiplicityPerp =.9;
			ApproximateEstimator.num_maxIterationsADMM = 1000;
			ApproximateEstimator.verboseLevel = 1;
			ApproximateEstimator.stepSizeADMM = .1;
			ApproximateEstimator.epsilon=.2;
			ApproximateEstimator.reg_separation = 10;
			ApproximateEstimator.str_criterion = 'approximate';
			% Construct an estimator for gossiping method
			GossipingEstimator = GossipingDecentralizedProjectionEstimator();
			GossipingEstimator.b_verbose=1;
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = 0;
			FastEstimator.reg_evalMultiplicityPerp =1;
			FastEstimator.reg_evalMultiplicityPar = .1;
			FastEstimator.reg_evalMultiplicityPerp =.9;
			FastEstimator.num_maxIterationsADMM = 1000;
			FastEstimator.epsilon=.2;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			% algorithm
			for ind_num=1:length(v_threshold)
				graphGen.threshold=v_threshold(ind_num);
				% Computing NMPE
				[m_Y(ind_num+counter,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges , FastEstimator);
				[ m_Y(ind_num+counter+1,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace, num_localExchanges , ApproximateEstimator);
				[ m_Y(ind_num+counter+2,:)] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
					dim_subspace , num_localExchanges , GossipingEstimator);
				counter=counter+2;
			end
			% Plot the results
			F = GFigure('m_X',1:num_localExchanges,'m_Y', m_Y,'ch_xlabel','Number of local exchanges',...
				'ch_ylabel','NMPE(H_l)','c_legend',{'Exact Projection,d_{max}=.25','Approximate Projection,d_{max}=.25','Gossiping,d_{max}=.25',...
				'Exact Projection, d_{max}=.35' ,'Approximate Projection,d_{max}=.35','Gossiping, d_{max}=.35'});
		end
		%% Figure 12
		% Monte Carlo simulation which represents the error between
		% designed filter and redesigned filter versus d_{max} to
		% show the effect of edge fluctuation
		% - Graph: Erdos-Renyi random graph
		function F = experiment_3012(obj,niter)
			% Set estimation problem parameters
			num_nodes =20;
			v_threshold=.35:0.05:0.8;
			v_dim_subspace =[2,3];
			num_localExchanges = num_nodes;
			counter=0;
			% Graph generator
			graphGen = WSNGraphGenerator();
			graphGen.num_nodes = num_nodes ;
			
			% Construct an estimator for fast decentralized method and set its parameters
			FastEstimator = FastDecentralizedProjectionEstimator();
			FastEstimator.str_criterion = 'approximate';
			FastEstimator.reg_energyPar = 0;
			FastEstimator.reg_energyPerp = 1;
			FastEstimator.reg_evalMultiplicityPar = 0.1;
			FastEstimator.reg_evalMultiplicityPerp = 0.9;
			FastEstimator.reg_separation=10;
			FastEstimator.num_maxIterationsADMM =1000;
			FastEstimator.verboseLevel = 1;
			FastEstimator.stepSizeADMM = .1;
			FastEstimator.epsilon=.1;
			
			num_failuresOne = 1;
			num_failuresTwo = 2;
			
			
			% simulation
			for ind_dim=1:length(v_dim_subspace)
				for ind_threshold=1:length(v_threshold)
					graphGen.threshold =v_threshold(ind_threshold);
					% Run the simulation
					[ NMPE ] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
						v_dim_subspace(ind_dim) , num_localExchanges , FastEstimator);
					[ NMPE_failuresOne] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
						v_dim_subspace(ind_dim) , num_localExchanges , FastEstimator, num_failuresOne);
					[ NMPE_failuresTwo] = FilterNMSEvsNumLocalExchanges( niter, graphGen ,...
						v_dim_subspace(ind_dim) , num_localExchanges , FastEstimator, num_failuresTwo);
					m_Y(counter+ind_dim,ind_threshold)= NMPE(end);
					m_Y(counter+ind_dim+1,ind_threshold)=NMPE_failuresOne(end);
					m_Y(counter+ind_dim+2,ind_threshold)=NMPE_failuresTwo(end);
				end
				counter=counter+2;
			end
			m_X=[v_threshold;v_threshold;v_threshold];
			m_X=[m_X;m_X];
			% Plot the results
			F = GFigure('m_X',m_X,'m_Y', m_Y,'ch_xlabel','p_{miss}',...
				'ch_ylabel','NMPE(H_l)','c_legend',...
				{'No node failure, r=2','One node failure, r=2','Two node failures, r=2','No node failure, r=3','One node failure, r=3','Two node failures, r=3'});
			
		end
		
		%% Figure 13
		% Num of distinct entries vs. number of blocks of the solutions to
		% a family of optimization problems designed to analyze the
		% relaxations in the paper. 
		function F = experiment_3013(obj,niter)
			function [m_A, v_b, v_sol] = block_mat_1(num_blocks, subblock_rows, subblock_cols)
				
				assert(subblock_rows+1 <= num_blocks*subblock_cols)
				m_A = [];
				for ind_block = 1:num_blocks
					block = randn(subblock_rows, subblock_cols);
					m_A = [m_A, block];
				end
				m_A = [m_A; ones(1, size(m_A,2))];
				v_b = [zeros(size(m_A,1)-1,1); 1];
				
				m_last_cols =  m_A(:,subblock_rows+1:end);
				m_reduced_A = [m_A(:,1:subblock_rows), sum(m_last_cols,2)];
				v_reduced_sol = m_reduced_A\v_b;
				num_last_cols = size(m_last_cols,2);
				v_sol = [v_reduced_sol(1:end-1);...
					v_reduced_sol(end)*ones(num_last_cols,1)];
				
				err = norm(m_A*v_sol - v_b);
				if err > 0.01
					err
					keyboard
				end
				
			end
			
			function [m_A, v_b, v_sol] = block_mat_2(num_blocks, subblock_rows, subblock_cols)
				
				%assert(subblock_size >= num_blocks + 1)
				assert(subblock_rows >= num_blocks)
				m_A = [];
				for ind_block = 1:num_blocks
					block = randn(subblock_rows, subblock_cols-1);
					block = [block, -sum(block,2)];
					below = zeros(num_blocks-1, size(block,2));
					if ind_block < num_blocks
						below(ind_block,1) = -1/ind_block;
					else
						below(:,1) = (1/ind_block); %* ones(num_blocks-1, size(block,2));
					end
					try
						block = [block;below];
					catch
						keyboard
					end
					
					m_A = [m_A, block];
				end
				%m_A = [m_A, [zeros(subblock_size,1); ones(size(below,1),1)]];
				m_A = [m_A; ones(1, size(m_A,2))];
				v_b = [zeros(size(m_A,1)-1,1); subblock_cols*sum(1:num_blocks)];
				
				% Sparse sol first
				v_subsol = m_A(:,1:subblock_cols)\v_b;
				v_sol = [v_subsol; zeros((num_blocks-1)*(subblock_cols),1)];
				%if (subblock_rows+2<subblock_cols) && (num_blocks > subblock_cols)
				%else
				if norm(m_A*v_sol - v_b) > 0.01
					v_sol = kron((1:num_blocks)',ones(subblock_cols,1));
				end
				
				err = norm(m_A*v_sol - v_b);
				if err > 0.01
					err
					keyboard
				end
				
			end
			
			function [m_A, v_b, v_sol] = block_mat_3(num_blocks, subblock_rows, subblock_cols)
				
				assert(subblock_rows >= num_blocks)
				%assert(subblock_cols >= num_blocks)
				m_A = [];
				for ind_block = 1:num_blocks
					block = randn(subblock_rows, subblock_cols-1);
					block = [block, -sum(block,2)];
					below = zeros(num_blocks-1, size(block,2));
					if ind_block < num_blocks
						below(ind_block,:) = -1/ind_block;
					else
						below(:) = (1/ind_block); %* ones(num_blocks-1, size(block,2));
					end
					try
						block = [block;below];
					catch
						keyboard
					end
					
					m_A = [m_A, block];
				end
				%m_A = [m_A, [zeros(subblock_size,1); ones(size(below,1),1)]];
				m_A = [m_A; ones(1, size(m_A,2))];
				v_b = [zeros(size(m_A,1)-1,1); subblock_cols*sum(1:num_blocks)];
				
				v_sol = kron((1:num_blocks)',ones(subblock_cols,1));
				
				err = norm(m_A*v_sol - v_b);
				if err > 0.01
					err
					keyboard
				end
				
			end
			
			function [num_opt, num_relaxed, num_max] = compare_num_distinct_elements(m_A, v_b, v_sol)
				v_relaxed_sol = DecentralizedProjectionExperiments.relaxMinDistinctEntrySol(m_A, v_b);
				num_relaxed = DecentralizedProjectionExperiments.UniqueTolerance(v_relaxed_sol,evalTol);
				num_opt = DecentralizedProjectionExperiments.UniqueTolerance(v_sol,evalTol);
				num_max = size(m_A,2);
			end
			
			evalTol = 1e-5;
			v_num_blocks = 3:10;
			
			function [num_opt, num_relaxed, num_max] = num_distinct_elements_vs_num_blocks(fun_mat)
				for ind=1:length(v_num_blocks)
					num_blocks = v_num_blocks(ind);
					%[m_A, v_b, v_sol] = block_mat_1(num_blocks, 6,6);
					%[m_A, v_b, v_sol] = block_mat_2(num_blocks, 2,10);
					%[m_A, v_b, v_sol] = block_mat_2(num_blocks, 8,8);
					%[m_A, v_b, v_sol] = fun_mat(num_blocks, 8,8);
					[m_A, v_b, v_sol] = fun_mat(num_blocks, 10,10);
					[num_opt(ind), num_relaxed(ind), num_max(ind)] = compare_num_distinct_elements(m_A, v_b, v_sol);
				end
			end
			
			[num_opt_1, num_relaxed_1, num_max_1] = num_distinct_elements_vs_num_blocks(@block_mat_1);
			[num_opt_2, num_relaxed_2, num_max_2] = num_distinct_elements_vs_num_blocks(@block_mat_2);
			[num_opt_3, num_relaxed_3, num_max_3] = num_distinct_elements_vs_num_blocks(@block_mat_3);
% 			
			F = GFigure('m_Y', [num_opt_1; num_relaxed_1; ...
				num_opt_2; num_relaxed_2; ...
				num_relaxed_3; num_max_3;],...
				'm_X',v_num_blocks,...
				'c_legend', {"Optimal D1", "Relaxed D1", ...
				"Optimal D2, D3", "Relaxed D2", ...
				"Relaxed D3", "Maximum D1, D2, D3"},...
				'c_styles', {"xr","-r","ob","-.b","--k",":g"},...
				'm_colorset', [1 0 0;1 0 0;0 .65 0; 0 .65 0;0 .65 0; 0 0 .9 ; .9 0 .9 ;.5 .5 0;0 .7 .7;...
				.5 0 1; 1 .5 0;  1 0 .5; 0 1 .5;0 1 0],...
				'ch_xlabel', "Number of blocks, B",...
				'ch_ylabel', "Number of distinct entries",...
				'ch_interpreter', "tex"...
				);
		end
	end
	
	
	methods(Static) % Methods to assist in the experiments
		
		function m_vander = maxFullColRankVanderMat(v_secondCol)
			% Vandermonde matrix with as many linearly indep. columns as
			% allowed by the finite precision arithmetic.
			%
			% Code can be optimized if this function is to be used
			% frequently.
			
			% Construction 1b: max number of columns with full rank
			num_evals = size(v_secondCol,1);
			m_vander = ones(num_evals,1);
			m_vander = m_vander/norm(m_vander);
			for ind_column = 2:num_evals
				v_oldColumn = m_vander(:,ind_column-1);
				v_newColumn = v_oldColumn.*v_secondCol;
				v_newColumn = v_newColumn/norm(v_newColumn);
				if rank([m_vander,v_newColumn]) == ind_column
					m_vander = [m_vander,v_newColumn];
				else
					break
				end
			end
			
			
		end
		
		
	end
	
end