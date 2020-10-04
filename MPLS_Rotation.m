%% Author: Yunpeng Shi
%% © Regents of the University of Minnesota. All rights reserved
%%------------------------------------------------
%% Message Passing Least Squares for Rotation Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Edge_Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). edge_num is the number of edges.
%% Cycle_Ind: 
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations
%% MPLS_parameters.max_iter: the number of iterations of MPLS
%% MPLS_parameters.beta: the sequence of reweighting parameter beta_t
%% MPLS_parameters.nsample: the number of sampled cycles per edge

%% Output:
%% SVec: Estimated corruption levels of all edges

%% Reference
%% [1] Yunpeng Shi and Gilad Lerman. "Message Passing Least Squares Framework and its Application to Rotation Synchronization" ICML 2020.




function[R_est, R_init] = MPLS_Rotation(Ind,RijMat,CEMP_parameters, MPLS_options)

    T=CEMP_parameters.max_iter;
    beta_vec=CEMP_parameters.reweighting_parameters;
    nsample = CEMP_parameters.nsample;
    T_beta = length(beta_vec);
    if T_beta<T
        beta_vec = [beta_vec,beta_vec(end)*(ones(1,T-T_beta))];
    end
    
    stop_threshold=MPLS_options.stop_threshold;
    maxIters = MPLS_options.max_iter;
    beta_mpls = MPLS_options.reweighting_parameters;
    T_beta_mpls = length(beta_mpls);
    if T_beta_mpls<maxIters
        beta_mpls = [beta_mpls,beta_mpls(end)*(ones(1,maxIters-T_beta_mpls))];
    end
    
    tau_mpls = MPLS_options.thresholding_parameters;
    T_tau_mpls = length(tau_mpls);
    if T_tau_mpls<maxIters
        tau_mpls = [tau_mpls,tau_mpls(end)*(ones(1,maxIters-T_tau_mpls))];
    end
    
    alpha_mpls = MPLS_options.cycle_info_ratio;
    T_alpha_mpls = length(alpha_mpls);
    if T_alpha_mpls<maxIters
        alpha_mpls = [alpha_mpls,alpha_mpls(end)*(ones(1,maxIters-T_alpha_mpls))];
    end
    
    
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');
    m=size(Ind_i,1);
    AdjMat = sparse(Ind_i,Ind_j,1,n,n);
    AdjMat = full(AdjMat + AdjMat');
    
    
    % start CEMP iterations as initialization
    
    disp('triangle sampling')
    %compute cycle inconsistency


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

    disp('Triangle Sampling Finished!')

   
    disp('Initializing')

    



    RijMat4d = zeros(3,3,n,n);
    % store relative rotations in 3x3xnxn tensor
    % construct edge index matrix (for 2d-to-1d index conversion)
    for l = 1:m
        i=Ind_i(l);j=Ind_j(l);
        RijMat4d(:,:,i,j)=RijMat(:,:,l);
        RijMat4d(:,:,j,i)=(RijMat(:,:,l))';
        IndMat(i,j)=l;
        IndMat(j,i)=-l;
    end
   


   

   
    Rki0 = zeros(3,3,m,nsample);
    Rjk0 = zeros(3,3,m,nsample);
    for l = IndPos
    Rki0(:,:,l,:) = RijMat4d(:,:,CoIndMat(:,l), Ind_i(l));
    Rjk0(:,:,l,:) = RijMat4d(:,:,Ind_j(l),CoIndMat(:,l));
    end
   

   
    
    Rki0Mat = reshape(Rki0,[3,3,m*nsample]);
    Rjk0Mat = reshape(Rjk0,[3,3,m*nsample]);
    Rij0Mat = reshape(kron(ones(1,nsample),reshape(RijMat,[3,3*m])), [3,3,m*nsample]);
   

    
    
    R_cycle0 = zeros(3,3,m*nsample);
    R_cycle = zeros(3,3,m*nsample);
    for j = 1:3
      R_cycle0 = R_cycle0 + bsxfun(@times,Rij0Mat(:,j,:),Rjk0Mat(j,:,:));
    end

    for j = 1:3
      R_cycle = R_cycle + bsxfun(@times,R_cycle0(:,j,:),Rki0Mat(j,:,:));
    end
    
    R_trace = (reshape(R_cycle(1,1,:)+R_cycle(2,2,:)+R_cycle(3,3,:), [m,nsample]))';
    S0Mat = abs(acos((R_trace-1)./2))/pi;

    

    SVec = mean(S0Mat,1);
    SVec(~IndPosbin)=1;
    disp('Initialization completed!')


    
    disp('Reweighting Procedure Started ...')


   for iter = 1:T
       
        % parameter controling the decay rate of reweighting function
        beta = beta_vec(iter);
        Ski = zeros(nsample, m);
        Sjk = zeros(nsample, m);
        for l = IndPos
            i = Ind_i(l); j=Ind_j(l);
            Ski(:,l) = SVec(abs(IndMat(i,CoIndMat(:,l))));
            Sjk(:,l) = SVec(abs(IndMat(j,CoIndMat(:,l))));
        end
        Smax = Ski+Sjk;
        % compute weight matrix (nsample by m)
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
        SMatij = sparse(Ind_j,Ind_i,SVec+1,n,n);
        [MST,~]=graphminspantree(SMatij);
        AdjTree = logical(MST+MST');
        
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
        
        % start MPLS iterations
        
        % transform the data format for the following Lie-Alegbraic
        % Averaging (LAA) solver --- matlab code written by Chatterjee
        RR = permute(RijMat, [2,1,3]); % relative rotations -- take transpose as the LAA code estimates R'
        Ind_T = Ind'; % indices matrix
        R_init = R_est; % use CEMP+MST as initialization for R
 
        
        % Formation of A matrix.
        Amatrix = Build_Amatrix(Ind_T);
        
        
        Q = R2Q(R_init); 
        QQ = R2Q(RR);
        
        

        score=inf;    Iteration=1;

        
        % initialization weights using (sij) estimated by CEMP
        Weights = (1./(SVec.^0.75)');

        
        
        weight_max = 1e4;
        weight_min = 1e-4;      
        thresh = quantile(SVec,tau_mpls(Iteration));
        Weights(Weights>weight_max)= weight_max;
        w=zeros(size(QQ,1),4);W=zeros(n,4);
        

        while((score>stop_threshold)&&(Iteration<maxIters))
             
  
        % one iteration of Weighted Lie-Algebraic Averaging                 
        [Q,W,B,score] = Weighted_LAA(Ind_T,Q,QQ,Amatrix,Weights);
 

            E=(Amatrix*W(2:end,2:4)-B); 
            ResVec = sqrt(sum(E.^2,2))/pi;

            Ski = zeros(nsample, m);
            Sjk = zeros(nsample, m);
            for l = IndPos
                i = Ind_i(l); j=Ind_j(l);
                Ski(:,l) = ResVec(abs(IndMat(i,CoIndMat(:,l))));
                Sjk(:,l) = ResVec(abs(IndMat(j,CoIndMat(:,l))));
            end
            Smax = Ski+Sjk;
            % compute weight matrix (nsample by m)
            WeightMat = exp(-beta_mpls(Iteration)*Smax);
            
            weightsum = sum(WeightMat,1);
            % normalize so that each column sum up to 1
            WeightMat = bsxfun(@rdivide,WeightMat,weightsum);
            HMat = WeightMat.*S0Mat;
            
            HVec = (sum(HMat,1))';
            RHVec = (1-alpha_mpls(Iteration))*ResVec + alpha_mpls(Iteration) * HVec;
           
            Weights = (1./(RHVec.^0.75));
            
           
            
            thresh = quantile(RHVec,tau_mpls(Iteration));
            Weights(Weights>weight_max)= weight_max;
            Weights(RHVec>thresh)=weight_min;
            
            fprintf('Iter %d: ||\x394R||= %f\n', Iteration, score); 
            Iteration = Iteration+1;
            

            
        end


        R_est=zeros(3,3,n);
        for i=1:size(Q,1)
            R_est(:,:,i)=q2R(Q(i,:));
        end
        
        
 
        
end

