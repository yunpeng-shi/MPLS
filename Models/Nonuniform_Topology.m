%% Author: Yunpeng Shi

%%------------------------------------------------
%% generation of the synthetic data with nonuniform topology
%%------------------------------------------------
%% Input Parameters: 
%% n: the number of the graph nodes
%% p: the probability of connecting a pair of vertices. G([n],E) is Erdos-Renyi graph G(n,p).
%% p_node_crpt: the probability of corrupting a node
%% p_node_edge: the probability of corrupting an edge
%% sigma_in: the noise level for inliers
%% sigma_out: the noise level for outliers
%% crpt_type (optional): choose 'uniform' or 'self-consistent', or 'adv'. The default choice is 'uniform'.

%% Output:
%% model_out.AdjMat: n by n adjacency matrix of the generated graph
%% model_out.Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). edge_num is the number of edges.
%% model_out.RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations
%% model_out.Rij_orig: 3 by 3 by edge_num tensor that stores the ground truth relative rotations
%% model_out.R_orig = R_orig: 3 by 3 by n tensor that stores the ground truth absolute rotations
%% model_out.ErrVec: the true corruption level of each edge
%% Reference
%% [1] Yunpeng Shi and Gilad Lerman. "Message Passing Least Squares Framework and its Application to Rotation Synchronization" ICML 2020.


function[model_out]=Nonuniform_Topology(n,p, p_node_crpt,p_edge_crpt, sigma_in, sigma_out, crpt_type)
    if ~exist('crpt_type','var')
        crpt_type = 'uniform';
    end
    
    G = rand(n,n) < p;
    G = tril(G,-1);
    % generate adjacency matrix
    AdjMat = G + G'; 
    [Ind_j, Ind_i] = find(G==1);
    Ind = [Ind_i,Ind_j];
    Ind_full = [Ind_j, Ind_i;Ind_i, Ind_j];
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
    IndMat = zeros(n,n);
    for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        Rij_orig(:,:,k)=R_orig(:,:,i)*(R_orig(:,:,j)');
        IndMat(i,j)=k;
        IndMat(j,i)=-k;
    end
    RijMat = Rij_orig;
    

    node_crpt = randperm(n);
    n_node_crpt = floor(n*p_node_crpt);
    node_crpt = node_crpt(1:n_node_crpt);
    crptInd = false(1,m);
    R_crpt = zeros(3,3,n); % only for self-consistent corruption

    for i = 1:n
        Q=randn(3);
        [U, ~, V]= svd(Q);
        S0 = diag([1,1,det(U*V')]);  
        R_crpt(:,:,i)=U*S0*V';
    end




    for i = node_crpt
        neighbor_cand = Ind_full(Ind_full(:,1)==i,2);
        neighbor_cand = reshape(neighbor_cand, 1, length(neighbor_cand));
        neighbor_crpt = randperm(length(neighbor_cand));
        n_neighbor = floor(p_edge_crpt * length(neighbor_cand));
        neighbor_crpt =    neighbor_crpt(1:n_neighbor);
        neighbor_crpt = neighbor_cand(neighbor_crpt);

        for j = neighbor_crpt 

            k = IndMat(i,j);
            crptInd(abs(k))=1;       
            Q=randn(3);
            [U, ~, V]= svd(Q);
            S0 = diag([1,1,det(U*V')]);  
            R0=U*S0*V';

            if strcmp(crpt_type,'uniform')
                if k>0
                    RijMat(:,:,k)= R0;    
                else
                    RijMat(:,:,-k)= R0';
                end

            end
            if strcmp(crpt_type,'self-consistent')
                if k>0
                RijMat(:,:,k)=R_crpt(:,:,i)*(R_crpt(:,:,j)');
                else
                RijMat(:,:,-k)=(R_crpt(:,:,i)*(R_crpt(:,:,j)'))';
                end
            end

            if strcmp(crpt_type,'adv')
                if k>0
                RijMat(:,:,k)=R_crpt(:,:,i)*R_orig(:,:,j)';
                else
                RijMat(:,:,-k)=(R_crpt(:,:,i)*R_orig(:,:,j)')';
                end
            end
        end     

    end
    
    
    noiseInd = ~crptInd;
    % indices of corrupted edges
    RijMat(:,:,noiseInd)= ...
    RijMat(:,:,noiseInd)+sigma_in*randn(3,3,sum(noiseInd)); % add noise

    RijMat(:,:,crptInd)= ...
    RijMat(:,:,crptInd)+sigma_out*randn(3,3,sum(crptInd)); % add noise
   



% project back to SO(3)    
    for k = 1:m
        [U, ~, V]= svd(RijMat(:,:,k));
        S0 = diag([1,1,det(U*V')]);
        RijMat(:,:,k) = U*S0*V';
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