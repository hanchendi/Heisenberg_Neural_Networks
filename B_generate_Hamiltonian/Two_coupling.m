function [H] = Two_coupling(p1,p2,sigma1,sigma2)

N=4;
I=eye(2);

if p1>1
    
    H=eye(2);
    for i=2:p1-1
        H=kron(H,I);
    end
    
    H=kron(H,sigma1);

    for i=p1+1:p2-1
        H=kron(H,I);
    end
    
    H=kron(H,sigma2);
    
    for i=p2+1:N
        H=kron(H,I);
    end
    
else
    
    H=sigma1;
    for i=2:p2-1
        H=kron(H,I);
    end
    
    H=kron(H,sigma2);
    
    for i=p2+1:N
        H=kron(H,I);
    end

end

