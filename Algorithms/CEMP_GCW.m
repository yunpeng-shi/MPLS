%% Author: Yunpeng Shi
%% 
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind

%% CEMP_parameters.max_iter: total # of iterations for CEMP
%% CEMP_parameters.reweighting: the sequence of reweighting parameter beta_t for CEMP
%% CEMP_parameters.nsample: # of cycles sampled per edge (also used in MPLS stage)



%% Output:
%% R_est: Estimated rotations (3x3xn)



function R_est = CEMP_GCW(Ind,RijMat,CEMP_parameters)

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
        SMat_sq = WeightMat.*S0Mat;
        % sij at current iteration
        SVec = sum(SMat_sq,1);
        SVec(~IndPosbin)=1;
        fprintf('Reweighting Iteration %d Completed!\n',iter)   
   end
        disp('Completed!')
    
    d=3;
    mat_size = ones(1,n)*d;
    cum_ind = [0,cumsum(mat_size)];
    Rij_blk = zeros(n*d);
    for k = 1:m
       i = Ind_i(k); j=Ind_j(k);
       Rij_blk((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1))= RijMat(:,:,k);    
    end

    Rij_blk = Rij_blk+Rij_blk';    
   
    SMat_sq = sparse(Ind_i, Ind_j, SVec, n, n);
    SMat_sq = SMat_sq + SMat_sq';
    Weights = exp(-beta.*SMat_sq).*AdjMat;
    Weights = diag(1./sum(Weights,2))*Weights;
    Weights = kron(Weights, ones(d));    
    RijW = Rij_blk.*Weights;
    clear 'Rij_blk';


    [V,~] = eigs(RijW,d,'la');
    V(:,1) = V(:,1)*sign(det(V(1:d,:))); % ensure det = 1
    R_est = zeros(d,d,n);
    for i=1:n
       Ri = V((cum_ind(i)+1):cum_ind(i+1), :); 
       [Ur,~,Vr] = svd(Ri);
       S0 = diag([ones(1,d-1),det(Ur*Vr')]);
       R_est(:,:,i) = Ur*S0*Vr';

    end  
        
        
end
