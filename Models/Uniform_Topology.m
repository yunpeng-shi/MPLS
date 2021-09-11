%% Author: Yunpeng Shi
%% © Regents of the University of Minnesota. All rights reserved
%%------------------------------------------------
%% generation of the synthetic data
%%------------------------------------------------
%% Input Parameters: 
%% n: the number of the graph nodes
%% p: the probability of connecting a pair of vertices. G([n],E) is Erdos-Renyi graph G(n,p).
%% q: the probability of corrupting an edge
%% sigma: the noise level (>0)
%% crpt_type (optional): choose 'uniform' or 'self-consistent'. The default choice is 'uniform'.

%% Output:
%% model_out.AdjMat: n by n adjacency matrix of the generated graph
%% model_out.Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). edge_num is the number of edges.
%% model_out.RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations
%% model_out.Rij_orig: 3 by 3 by edge_num tensor that stores the ground truth relative rotations
%% model_out.R_orig = R_orig: 3 by 3 by n tensor that stores the ground truth absolute rotations
%% model_out.ErrVec: the true corruption level of each edge
%% Reference
%% [1] Yunpeng Shi and Gilad Lerman. "Message Passing Least Squares Framework and its Application to Rotation Synchronization" ICML 2020.


function[model_out]=Uniform_Topology(n,p,q,sigma,crpt_type)
    if ~exist('crpt_type','var')
        crpt_type = 'uniform';
    end
    
    G = rand(n,n) < p;
    G = tril(G,-1);
    % generate adjacency matrix
    AdjMat = G + G'; 
    [Ind_j, Ind_i] = find(G==1);
    Ind = [Ind_i,Ind_j];
    m = length(Ind_i);

    %generate rotation matrices
    R_orig = zeros(3,3,n);

    for i = 1:n
        Q=randn(3);
        [U, ~, V]= svd(Q);
        S0 = diag([1,1,det(U*V')]);  
        R_orig(:,:,i)=U*S0*V';
    end

    Rij_orig = zeros(3,3,m);
    for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        Rij_orig(:,:,k)=R_orig(:,:,i)*(R_orig(:,:,j)');
    end
    RijMat = Rij_orig;
    noiseIndLog = rand(1,m)>=q;
    % indices of corrupted edges
    corrIndLog = logical(1-noiseIndLog);
    noiseInd=find(noiseIndLog);
    corrInd=find(corrIndLog);
    RijMat(:,:,noiseInd)= ...
    RijMat(:,:,noiseInd)+sigma*randn(3,3,length(noiseInd)); % add noise
    % project back to SO(3)
    for k = noiseInd
        [U, ~, V]= svd(RijMat(:,:,k));
        S0 = diag([1,1,det(U*V')]);
        RijMat(:,:,k) = U*S0*V';
    end    
    
    R_corr = zeros(3,3,n); % only for self-consistent corruption

    for i = 1:n
        Q=randn(3);
        [U, ~, V]= svd(Q);
        S0 = diag([1,1,det(U*V')]);  
        R_corr(:,:,i)=U*S0*V';
    end

    if strcmp(crpt_type,'uniform')
        for k = corrInd
            Q=randn(3);
            [U, ~, V]= svd(Q);
            S0 = diag([1,1,det(U*V')]);
            RijMat(:,:,k) = U*S0*V';    % corrupted relative rotations that follows uniform distribtion
        end
    else
        for k = corrInd
            i=Ind_i(k); j=Ind_j(k); 
            Q=R_corr(:,:,i)*(R_corr(:,:,j)')+sigma*randn(3,3);
            [U, ~, V]= svd(Q);
            S0 = diag([1,1,det(U*V')]);
            RijMat(:,:,k) = U*S0*V'; % corrupted relative rotations that are self-conisistent
        end
    end
        

    R_err = zeros(3,3,m);
    for j = 1:3
      R_err = R_err + bsxfun(@times,Rij_orig(:,j,:),RijMat(:,j,:));
    end


    R_err_trace = (reshape(R_err(1,1,:)+R_err(2,2,:)+R_err(3,3,:), [m,1]))';
    ErrVec = abs(acos((R_err_trace-1)./2))/pi; % computing geodesic distance from the ground truth
    
    
    model_out.AdjMat = AdjMat;
    model_out.Ind = Ind;
    model_out.RijMat = RijMat;
    model_out.Rij_orig = Rij_orig;
    model_out.R_orig = R_orig;
    model_out.ErrVec = ErrVec;
    
    
end
