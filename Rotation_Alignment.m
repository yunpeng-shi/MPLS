%% Rotation Alignment
%% This function finds rotation R that minimizes sum_i ||R_i^* - R_i*R||_F^2
%% Input: 
%% Ri  : d-by-d-n matrix of rotations to be registered ({Ri} above) 
%% Ri0 : d-by-d-n matrix of reference rotations ({Ri^0} above) 
%% 
%% Output:
%% R_fit   : Ri registered to Ri0 (R_fit = R_opt*Ri)
%% MSE_rot : Mean squared error 
%% R_opt   : d-by-d Optimal point of (1)
%%*************************************************************************

function [R_fit, R_opt, mean_error, median_error] = Rotation_Alignment(Ri, Ri0)

d = size(Ri0,1); n = size(Ri0,3);
Q = zeros(d,d);

for k = 1:n
    A_est_i = Ri(:,:,k);
    A_orig_i = Ri0(:,:,k);
    Q = Q + A_est_i'*A_orig_i;
end

[Uq,~,Vq] = svd(Q);
R_opt = Uq*[eye(d-1) zeros(d-1,1); zeros(1,d-1) det(Uq*Vq')]*Vq';


MSEVec = zeros(1,n);
R_fit = zeros(d,d,n);
for k = 1:n
    R_fit(:,:,k) = Ri(:,:,k)*R_opt;
    %MSE_rot = MSE_rot + norm(Ri0(:,:,k) - R_fit(:,:,k), 'fro')^2;
    R_tr = trace(Ri0(:,:,k)*(R_fit(:,:,k))');
    MSEVec(k) =  abs(acos((R_tr-1)./2))/pi*180;
end
mean_error = mean(MSEVec);
median_error = median(MSEVec);

return
