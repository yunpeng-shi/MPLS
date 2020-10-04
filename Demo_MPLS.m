
% generate synthetic data with uniform corruption

model_out = Rotation_Graph_Generation(200,0.5,0.85,0,'uniform');

Ind = model_out.Ind;
RijMat = model_out.RijMat;
ErrVec = model_out.ErrVec;
R_orig = model_out.R_orig;

% set CEMP defult parameters
CEMP_parameters.max_iter = 6;
CEMP_parameters.reweighting_parameters = 2.^((1:6)-1);
CEMP_parameters.nsample = 50;

% set MPLS default parameters and options
MPLS_options.stop_threshold = 1e-3;
MPLS_options.max_iter = 100;
MPLS_options.reweighting_parameters = CEMP_parameters.reweighting_parameters(end);
MPLS_options.thresholding_parameters = [0.95,0.9,0.85,0.8];
MPLS_options.cycle_info_ratio = 1./((1:MPLS_options.max_iter)+1);


% run MPLS
[R_MPLS, R_CEMP_MST] = MPLS_Rotation(Ind,RijMat,CEMP_parameters, MPLS_options);

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

