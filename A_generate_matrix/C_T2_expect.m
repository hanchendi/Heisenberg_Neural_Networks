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

T2=zeros(N^2,N^2,2*N,2*N);
for i1=1:N^2
    for i2=1:N^2
        for i3=1:2*N
            for i4=1:2*N
            
                x1=X1(:,:,i3);
                x1=reshape(x1,[1,N]);
                x2=X2(:,:,i1);
                x2=reshape(x2,[N,N]);
                x3=X2(:,:,i2);
                x3=reshape(x3,[N,N]);
                x4=X3(:,i4,:);
                x4=reshape(x4,[N,1]);
            
                T2(i1,i2,i3,i4)=x1*sqrt(-1)*(x2*x3-x3*x2)*x4;
            
            end
        end
    end
    disp([1,i1/N^2])
end

for i1=1:N^2
    for i2=1:N^2
        for i3=1:2*N
            for i4=1:i3
                if T2(i1,i2,i3,i4)==-T2(i1,i2,i4,i3) && T2(i1,i2,i3,i4)~=0
                    T2(i1,i2,i3,i4)=0;
                    T2(i1,i2,i4,i3)=0;
                end
            end
        end
    end
    disp([2,i1/N^2])
end

save('T2.mat','T2')