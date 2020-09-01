rng(3);
model_out = rotation_graph_generation(200,0.5,0.5,0,'uniform');

Ind = model_out.Ind;
RijMat = model_out.RijMat;
ErrVec = model_out.ErrVec;

CEMP_parameters.max_iter = 6;
CEMP_parameters.beta = 2.^((1:6)-1);
CEMP_parameters.nsample = 50;

SVec = CEMP_rotation(Ind,RijMat,CEMP_parameters);
plot(ErrVec,SVec,'b.');
title('Scatter Plot of s_{ij}^* v.s. s_{ij,T}')
xlabel('s_{ij}^*') 
ylabel('s_{ij,T}') 