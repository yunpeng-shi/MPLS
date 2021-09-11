%% Rotation Alignment
%% This function finds rotation R that minimizes sum_i ||R_i^* - R_i*R||_F^2
%% Input: 
%% R_est  : d-by-d-n matrix of rotations to be registered) 
%% R_gt : d-by-d-n matrix of reference rotations) 
%% 
%% Output:
%% R_out   : R_est registered to R_gt using R_out = R_est*R_align
%% Mean and Median error (geodesic distance in degree) of registration
%% R_align: The right rotation for alignment (the only variable to be optimized)


function [R_out, R_align, mean_error, median_error] = Rotation_Alignment(R_est, R_gt)

    d = size(R_gt,1); n = size(R_gt,3);
    A = zeros(d,d);

    for k = 1:n
        U = R_est(:,:,k);
        V = R_gt(:,:,k);
        A = A + U'*V;
    end

    [U1,~,V1] = svd(A);
    R_align = U1*[eye(d-1) zeros(d-1,1); zeros(1,d-1) det(U1*V1')]*V1';


    MSEVec = zeros(1,n);
    R_out = zeros(d,d,n);
    for k = 1:n
        R_out(:,:,k) = R_est(:,:,k)*R_align;
        R_tr = trace(R_gt(:,:,k)*(R_out(:,:,k))');
        MSEVec(k) =  abs(acos((R_tr-1)./2))/pi*180;
    end
    mean_error = mean(MSEVec);
    median_error = median(MSEVec);

end
