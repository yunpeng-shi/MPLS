%% Author: Yunpeng Shi
%% © Regents of the University of Minnesota. All rights reserved
%%------------------------------------------------
%% Cycle-Edge Message Passing for Rotation Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations
%% CEMP_parameters.max_iter: the number of iterations of CEMP
%% CEMP_parameters.beta: the sequence of reweighting parameter beta_t
%% CEMP_parameters.nsample: the number of sampled cycles per edge

%% Output:
%% SVec: Estimated corruption levels of all edges

%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing" arXiv preprint, 2019
%% [2] Yunpeng Shi and Gilad Lerman. "Message Passing Least Squares Framework and its Application to Rotation Synchronization" ICML 2020.




function[SVec] = CEMP_rotation(Ind,RijMat,CEMP_parameters)
    T=CEMP_parameters.max_iter;
    beta_vec=CEMP_parameters.reweighting_parameters;
    nsample = CEMP_parameters.nsample;
    T_beta = length(beta_vec);
    if T_beta<T
        beta_vec = [beta_vec,beta_vec(end)*(ones(1,T-T_beta))];
    end
            
        Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');
    m=size(Ind_i,1);
    AdjMat = sparse(Ind_i,Ind_j,1,n,n);
    AdjMat = full(AdjMat + AdjMat');
    
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
        
       
  

end