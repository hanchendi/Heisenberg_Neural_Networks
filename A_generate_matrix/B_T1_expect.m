clear
clc

N=2^4;
X1=zeros(1,N,2*N);
X3=zeros(N,2*N,1);
for i=1:N
    X1(1,i,2*i-1)=1;
    X1(1,i,2*i)=-sqrt(-1);
    X3(i,2*i-1,1)=1;
    X3(i,2*i,1)=sqrt(-1);
end

X2=zeros(N,N,N^2);

o=1;
for i=1:N
    for j=i:N
        
        if i==j
            X2(i,j,o)=1;
            o=o+1;
        else
            X2(i,j,o)=1;
            X2(i,j,o+1)=sqrt(-1);
            o=o+2;
        end
        
    end
end

for i=1:N
    for j=1:i-1
        for k=1:N^2
            X2(i,j,k)=conj(X2(j,i,k));
        end
    end
end

T1=zeros(N^2,2*N,2*N);
for i=1:N^2
    for j=1:2*N
        for k=1:2*N
            
            x1=X1(:,:,j);
            x1=reshape(x1,[1,N]);
            x2=X2(:,:,i);
            x2=reshape(x2,[N,N]);
            x3=X3(:,k,:);
            x3=reshape(x3,[N,1]);
            
            T1(i,j,k)=x1*x2*x3;
            
        end
    end
    disp([1,i/N^2])
end

for i=1:N^2
    for j=1:2*N
        for k=1:2*N
            if T1(i,j,k)==-T1(i,k,j) && T1(i,j,k)~=0
                T1(i,j,k)=0;
                T1(i,k,j)=0;
            end
        end
    end
    disp([2,i/N^2])
end

save('T1.mat','T1')