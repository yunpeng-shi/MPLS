%% Original Author: AVISHEK CHATTERJEE
%% Taken from the code of AVISHEK CHATTERJEE for Lie-Algebraic Averaging
%% This code aims to convert blocks of rotations to blocks of quaternions.
%% Reference: "Efficient and Robust Large-Scale Rotation Averaging." by Avishek Chatterjee, Venu Madhav Govindu.


function[Qtn]=R2Q(Rot)
        
        Qtn=[Rot(1,1,:)+Rot(2,2,:)+Rot(3,3,:)-1, Rot(3,2,:)-Rot(2,3,:),Rot(1,3,:)-Rot(3,1,:),Rot(2,1,:)-Rot(1,2,:)]/2;
        Qtn=reshape(Qtn,4,size(Qtn,3),1)';
        Qtn(:,1)=sqrt((Qtn(:,1)+1)/2);
        Qtn(:,2:4)=(Qtn(:,2:4)./repmat(Qtn(:,1),[1,3]))/2;

       
end