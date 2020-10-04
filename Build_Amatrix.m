function[Amatrix] = Build_Amatrix(I)
        m=size(I,2);
        N=max(max(I));
        i=[[1:m];[1:m]];i=i(:);
        j=I(:);
        s=repmat([-1;1],[m,1]);
        k=(j~=1);
        Amatrix=sparse(i(k),j(k)-1,s(k),m,N-1);
end