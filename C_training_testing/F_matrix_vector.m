function [M_V] = F_matrix_vector(M)

    M_V=zeros(256,1);
	o=1;
    for i1=1:16
        for i2=i1:16
            if i1==i2
                M_V(o)=real(M(i1,i2));
                o=o+1;
            else
                M_V(o)=(real(M(i1,i2)+M(i2,i1)))/2;
                o=o+1;
                M_V(o)=(imag(M(i1,i2)-M(i2,i1)))/2;
                o=o+1;
            end
        end
    end
    
end

