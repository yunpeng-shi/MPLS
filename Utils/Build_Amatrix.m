%% Original Author: AVISHEK CHATTERJEE
%% Taken from the code of AVISHEK CHATTERJEE for Lie-Algebraic Averaging
%% This code aims to build the A matrix for the least squares solver min||Ax-b||^2 for rotation averaging
%% Reference: "Efficient and Robust Large-Scale Rotation Averaging." by Avishek Chatterjee, Venu Madhav Govindu.

function[Amatrix] = Build_Amatrix(I)
        m=size(I,2);
        N=max(max(I));
        i=[[1:m];[1:m]];i=i(:);
        j=I(:);
        s=repmat([-1;1],[m,1]);
        k=(j~=1);
        Amatrix=sparse(i(k),j(k)-1,s(k),m,N-1);
end