function [H] = Self_coupling(p,sigma)

N=4;
I=eye(2);

if p>1
    
    H=eye(2); 
    for i=2:p-1
        H=kron(H,I);
    end

    H=kron(H,sigma);

    for i=p+1:N
        H=kron(H,I);
    end
    
else
    
    H=sigma;
    for i=2:N
        H=kron(H,I);
    end

end

