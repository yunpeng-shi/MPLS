rng(3)
% generate synthetic data with uniform corruption

model_out = Rotation_Graph_Generation(200,0.5,0.85,0,'uniform');

Ind = model_out.Ind; % matrix of edge indices (m by 2)
RijMat = model_out.RijMat; % given corrupted and noisy relative rotations
ErrVec = model_out.ErrVec; % ground truth corruption levels
R_orig = model_out.R_orig; % ground truth rotations

% set CEMP defult parameters
CEMP_parameters.max_iter = 6;
CEMP_parameters.reweighting = 2.^((1:6)-1);
CEMP_parameters.nsample = 50;

% set MPLS default parameters and options
MPLS_parameters.stop_threshold = 1e-3;
MPLS_parameters.max_iter = 100;
MPLS_parameters.reweighting = CEMP_parameters.reweighting(end);
MPLS_parameters.thresholding = [0.95,0.9,0.85,0.8];
MPLS_parameters.cycle_info_ratio = 1./((1:MPLS_parameters.max_iter)+1);

% For dense graphs with sufficient uncorrupted 3-cycles for all edges, 
% the following parameters may work even better: reweighting parameter can gradually
% increase (in ICML paper we fix beta=32 for MPLS). One can increasingly weigh 3-cycle
% consistency information and ignore residual nformation (in ICML paper we
% gradually ignore cycle information and weigh residual more). 

% MPLS_parameters.reweighting = 0.1*1.5.^((1:15)-1);
% MPLS_parameters.cycle_info_ratio = 1-1./((1:MPLS_parameters.max_iter)+1);



% run MPLS
[R_MPLS, R_CEMP_MST] = MPLS_Rotation(Ind,RijMat,CEMP_parameters, MPLS_parameters);

% rotation alignment for evaluation
[~, ~, mean_error_CEMP_MST, median_error_CEMP_MST] = Rotation_Alignment(R_CEMP_MST, R_orig);
[~, ~, mean_error_MPLS, median_error_MPLS] = Rotation_Alignment(R_MPLS, R_orig);



% Report estimation error
sz = [2 3];
varTypes = {'string','double','double'};
varNames = {'Algorithms','MeanError','MedianError'};
Results = table('Size',sz,'VariableTypes',varTypes, 'VariableNames',varNames);
Results(1,:)={'CEMP+MST', mean_error_CEMP_MST, median_error_CEMP_MST};
Results(2,:)={'MPLS', mean_error_MPLS, median_error_MPLS}

