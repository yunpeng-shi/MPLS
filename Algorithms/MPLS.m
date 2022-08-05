%% Author: Yunpeng Shi
%% © Regents of the University of Minnesota. All rights reserved
%%------------------------------------------------
%% Message Passing Least Squares for Rotation Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind

%% CEMP_parameters.max_iter: total # of iterations for CEMP
%% CEMP_parameters.reweighting: the sequence of reweighting parameter beta_t for CEMP
%% CEMP_parameters.nsample: # of cycles sampled per edge (also used in MPLS stage)

%% MPLS_parameters.stop_threshold: stopping criterion
%% MPLS_parameters.max_iter: the maximal number of iterations of MPLS
%% MPLS_parameters.reweighting: the sequence of reweighting parameter beta_t for MPLS
%% MPLS_parameters.thresholding: the sequence of thresholding parameters tau_t
%% MPLS_parameters.cycle_info_ratio: the coefficient alpha_t of cycle-consistency information.



%% Output:
%% R_est: Estimated rotations (3x3xn)
%% R_init: Initialized rotations by CEMP+MST (3x3xn)

%% Reference
%% [1] Yunpeng Shi and Gilad Lerman. "Message Passing Least Squares Framework and its Application to Rotation Synchronization" ICML 2020.


function[R_est, R_init] = MPLS(Ind,RijMat,CEMP_parameters, MPLS_parameters)

    %CEMP parameters
    T=CEMP_parameters.max_iter; 
    beta_cemp=CEMP_parameters.reweighting;
    nsample = CEMP_parameters.nsample;
    T_beta = length(beta_cemp);
    if T_beta<T
        % if the reweighting parameter vector is short, then the rest of
        % missing elements are set to constant
        beta_cemp = [beta_cemp,beta_cemp(end)*(ones(1,T-T_beta))]; 
    end
    
    % MPLS paramters
    stop_threshold=MPLS_parameters.stop_threshold;
    maxIters = MPLS_parameters.max_iter;
    beta_mpls = MPLS_parameters.reweighting;
    T_beta_mpls = length(beta_mpls);
    if T_beta_mpls<maxIters
        beta_mpls = [beta_mpls,beta_mpls(end)*(ones(1,maxIters-T_beta_mpls))];
    end
    
    tau_mpls = MPLS_parameters.thresholding;
    T_tau_mpls = length(tau_mpls);
    if T_tau_mpls<maxIters
        tau_mpls = [tau_mpls,tau_mpls(end)*(ones(1,maxIters-T_tau_mpls))];
    end
    
    alpha_mpls = MPLS_parameters.cycle_info_ratio;
    T_alpha_mpls = length(alpha_mpls);
    if T_alpha_mpls<maxIters
        alpha_mpls = [alpha_mpls,alpha_mpls(end)*(ones(1,maxIters-T_alpha_mpls))];
    end
    
    % building the graph   
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');
    m=size(Ind_i,1);
    AdjMat = sparse(Ind_i,Ind_j,1,n,n); % Adjacency matrix
    AdjMat = full(AdjMat + AdjMat');
        
    % start CEMP iterations as initialization   
    disp('sampling 3-cycles')
    % Matrix of codegree:
    % CoDeg(i,j) = 0 if i and j are not connected, otherwise,
    % CoDeg(i,j) = # of vertices that are connected to both i and j
    CoDeg = (AdjMat*AdjMat).*AdjMat;
    AdjPos = AdjMat;
    % label positive codegree elements as -1
    AdjPos(CoDeg>0)=-1;
    AdjPosLow = tril(AdjPos);
    AdjPosLow(AdjPosLow==0)=[];
    % find 1d indices of edges with positive codegree
    IndPos = find(AdjPosLow<0);
    IndPosbin = zeros(1,m);
    IndPosbin(IndPos)=1;
    % CoIndMat(:,l)= triangles sampled that contains l-th edge
    % e.g. l=(3,5), then CoIndMat(:,l)=[2,9,8,...] means that...
    % triangles 352, 359, 358,... are sampled
    for l = IndPos
        i = Ind_i(l); j = Ind_j(l);
       CoIndMat(:,l)= datasample(find(AdjMat(:,i).*AdjMat(:,j)), nsample);
    end

    disp('Sampling Finished!')
    disp('Initializing')

    RijMat4d = zeros(3,3,n,n);
    for l = 1:m
        i=Ind_i(l);j=Ind_j(l);
        RijMat4d(:,:,i,j)=RijMat(:,:,l); % store relative rotations in 3x3xnxn tensor
        RijMat4d(:,:,j,i)=(RijMat(:,:,l))';
        IndMat(i,j)=l; % construct edge index matrix (for 2d-to-1d index conversion)
        IndMat(j,i)=-l;
    end
   
    % start computing cycle-inconsistency
    Rki0 = zeros(3,3,m,nsample); % Rki0(:,:,l,s) is Rki if the s-th sampled cycle of edge ij (whose 1-d index is l) is ijk
    Rjk0 = zeros(3,3,m,nsample);
    for l = IndPos
    Rki0(:,:,l,:) = RijMat4d(:,:,CoIndMat(:,l), Ind_i(l));
    Rjk0(:,:,l,:) = RijMat4d(:,:,Ind_j(l),CoIndMat(:,l));
    end
 
    % reshape above matrices for easier multiplication
    Rki0Mat = reshape(Rki0,[3,3,m*nsample]);
    Rjk0Mat = reshape(Rjk0,[3,3,m*nsample]);
    Rij0Mat = reshape(kron(ones(1,nsample),reshape(RijMat,[3,3*m])), [3,3,m*nsample]);
    
    R_cycle0 = zeros(3,3,m*nsample);
    R_cycle = zeros(3,3,m*nsample);
    for j = 1:3
      R_cycle0 = R_cycle0 + bsxfun(@times,Rij0Mat(:,j,:),Rjk0Mat(j,:,:));
    end
    for j = 1:3
      R_cycle = R_cycle + bsxfun(@times,R_cycle0(:,j,:),Rki0Mat(j,:,:));  % R_cycle(:,:,s) stores Rij*Rjk*Rki for that cycle
    end 
    R_trace = (reshape(R_cycle(1,1,:)+R_cycle(2,2,:)+R_cycle(3,3,:), [m,nsample]))'; % R_trace(l,s) stores Tr(Rij*Rjk*Rki) for that cycle
    S0Mat = abs(acos((R_trace-1)./2))/pi;   % S0Mat(l,s) stores d(Rij*Rjk*Rki, I) for that cycle. d is the normalized geodesic distance.
    SVec = mean(S0Mat,1); % initialized corruption level estimates s_{ij,0}
    SVec(~IndPosbin)=1; % If there is no 3-cycle for that edge, then set its sij as 1 (the largest possible).
    disp('Initialization completed!')  
    disp('Reweighting Procedure Started ...')

   for iter = 1:T     
        % parameter controling the decay rate of reweighting function
        beta = beta_cemp(iter);
        Ski = zeros(nsample, m);
        Sjk = zeros(nsample, m);
        for l = IndPos
            i = Ind_i(l); j=Ind_j(l);
            Ski(:,l) = SVec(abs(IndMat(i,CoIndMat(:,l))));
            Sjk(:,l) = SVec(abs(IndMat(j,CoIndMat(:,l))));
        end
        Smax = Ski+Sjk;
        % compute cycle-weight matrix (nsample by m)
        WeightMat = exp(-beta*Smax);
        weightsum = sum(WeightMat,1);
        % normalize so that each column sum up to 1
        WeightMat = bsxfun(@rdivide,WeightMat,weightsum);
        SMat = WeightMat.*S0Mat;
        % sij at current iteration
        SVec = sum(SMat,1);
        SVec(~IndPosbin)=1;
        fprintf('Reweighting Iteration %d Completed!\n',iter)   
   end
        disp('Completed!')
        
        % find the minimum spanning tree given SVec-weighted graph
        disp('Building minimum spanning tree ...')
        SMatij = sparse(Ind_j,Ind_i,SVec+1,n,n);
%         [MST,~]=graphminspantree(SMatij);  This is the legacy code for matlab versions before 2020a
%         AdjTree = logical(MST+MST');
        
        G = graph(SMatij,'lower'); % These 3 lines are compatible with new matlab versions
        Tree = minspantree(G);
        AdjTree = adjacency(Tree);
        
        % compute Ri by multiplying Rij along the spanning tree
        rootnodes = 1;
        added=zeros(1,n);
        R_est = zeros(3,3,n);
        R_est(:,:,rootnodes)=eye(3);
        added(rootnodes)=1;
        newroots = [];
                
        while sum(added)<n
            for node_root = rootnodes
                leaves = find((AdjTree(node_root,:).*(1-added))==1);
                newroots = [newroots, leaves];
                for node_leaf=leaves
                    edge_leaf = IndMat(node_leaf,node_root);
                    if edge_leaf>0
                        R_est(:,:,node_leaf)=RijMat(:,:,abs(edge_leaf))*R_est(:,:,node_root);
                    else
                        R_est(:,:,node_leaf)=(RijMat(:,:,abs(edge_leaf)))'*R_est(:,:,node_root);
                    end
                    added(node_leaf)=1;
                end
            end
            rootnodes = newroots;
        end
        
        
        % start MPLS procedure
        
        % transform the data format for the following Lie-Alegbraic
        % Averaging (LAA) solver that was written by AVISHEK CHATTERJEE
        RR = permute(RijMat, [2,1,3]); % relative rotations -- take transpose as the original LAA code estimates R' in our setting
        Ind_T = Ind'; % indices matrix
        R_init = R_est; % use CEMP+MST as initialization for R        
        % Formation of A matrix.
        Amatrix = Build_Amatrix(Ind_T); % A matrix for least squares solver (by AVISHEK CHATTERJEE)      
        Q = R2Q(R_init); % Transfer to quoternion representation (by AVISHEK CHATTERJEE)  
        QQ = R2Q(RR); % Transfer to quoternion representation (by AVISHEK CHATTERJEE)  
        score=inf;    Iteration=1;
     
        % initialization edge weights using (sij) estimated by CEMP
        Weights = (1./(SVec.^0.75)');
        weight_max = 1e4; % weights cannot be too large nor too small (for numerical stability and graph connectivity)
        weight_min = 1e-4;      
        Weights(Weights>weight_max)= weight_max;
       
        disp('Rotation Initialized!')
        disp('Start MPLS reweighting ...')
        % start MPLS reweighting iterations
        while((score>stop_threshold)&&(Iteration<maxIters))
             % one iteration of Weighted Lie-Algebraic Averaging (by AVISHEK CHATTERJEE)                
            [Q,W,B,score] = Weighted_LAA(Ind_T,Q,QQ,Amatrix,Weights); 
            E=(Amatrix*W(2:end,2:4)-B); 
            ResVec = sqrt(sum(E.^2,2))/pi; % normalized residual rij for all edges
            Ski = zeros(nsample, m);
            Sjk = zeros(nsample, m);
            for l = IndPos
                i = Ind_i(l); j=Ind_j(l);
                Ski(:,l) = ResVec(abs(IndMat(i,CoIndMat(:,l))));
                Sjk(:,l) = ResVec(abs(IndMat(j,CoIndMat(:,l))));
            end
            Smax = Ski+Sjk;
            % compute cycle-weight matrix (nsample by m)
            WeightMat = exp(-beta_mpls(Iteration)*Smax);        
            weightsum = sum(WeightMat,1);
            % normalize so that each column sum up to 1
            WeightMat = bsxfun(@rdivide,WeightMat,weightsum);
            HMat = WeightMat.*S0Mat;
            HVec = (sum(HMat,1))'; % compute hij for each edge
            % The following commented line is optional:  
            % HVec(~IndPosbin)=1;
            RHVec = (1-alpha_mpls(Iteration))*ResVec + alpha_mpls(Iteration) * HVec; % convex combination of rij and hij     
            Weights = (1./(RHVec.^0.75)); % compute edge weights
            % additional truncation for edge weights
            thresh = quantile(RHVec,tau_mpls(Iteration));
            Weights(Weights>weight_max)= weight_max;
            Weights(RHVec>thresh)=weight_min;
            % report the change of estimated rotations (stop when the change is small) 
            fprintf('Iter %d: ||\x394R||= %f\n', Iteration, score); 
            Iteration = Iteration+1;        
        end
        % transform from quaternion and return the estimated rotations
        R_est=zeros(3,3,n);
        for i=1:size(Q,1)
            R_est(:,:,i)=q2R(Q(i,:));
        end 
        disp('DONE!')
        
end
