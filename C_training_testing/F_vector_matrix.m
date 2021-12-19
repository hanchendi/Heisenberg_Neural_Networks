function [M] = F_vector_matrix(M_V)

    M=zeros(16,16);
	o=1;
    for i1=1:16
        for i2=i1:16
            if i1==i2
                M(i1,i2)=M_V(o);
                o=o+1;
            else
                M(i1,i2)=M_V(o)+M_V(o+1)*sqrt(-1);
                o=o+2;
            end
        end
    end
    
    for i1=1:16
        for i2=1:i1-1
            M(i1,i2)=conj(M(i2,i1));
        end
    end
end

