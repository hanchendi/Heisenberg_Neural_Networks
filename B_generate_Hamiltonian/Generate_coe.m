clear
clc

% Please see Eq.(9)
% c_self: self coupling coefficient
% The first index represent h^{(1,2)}
% The second term represent 4 spins
% The third term represent sigma_{x,y,z}

% c_coupling: two body interaction term
% c_time: Eq.(10)

for s=1:1
    c_self=rand(2,4,3);
    c_coupling=rand(2,4,9);
    c_time=rand(1,3);
    
    w=zeros(4,2);
    w(1,:)=[1 2];
    w(2,:)=[2 3];
    w(3,:)=[3 4];
    w(4,:)=[1 4];

    sigma_gather=zeros(3,2,2);
    sigma_gather(1,:,:)=[0 1; 1 0];
    sigma_gather(2,:,:)=[0 -sqrt(-1); sqrt(-1) 0];
    sigma_gather(3,:,:)=[1 0; 0 -1];

    H_Matrix=zeros(2,16,16);
    for i1=1:2
    
        %%%%%%%%%%%%%%%
        % Self coupling
        %%%%%%%%%%%%%%%
    
        H=zeros(16,16);
        for i2=1:4   % number of spin
            for i3=1:3 % xyz
                sigma=sigma_gather(i3,:,:);
                sigma=reshape(sigma,[2,2]);
                H=H+c_self(i1,i2,i3)*Self_coupling(i2,sigma);
            end
        end
    
        %%%%%%%%%%%%%%%
        % Two coupling
        %%%%%%%%%%%%%%%    
    
        for i2=1:4
        
            o=1;
            for j1=1:3 % xyz
                for j2=1:3 % xyz
                
                    sigma1=sigma_gather(j1,:,:);
                    sigma1=reshape(sigma1,[2,2]);
                    sigma2=sigma_gather(j2,:,:);
                    sigma2=reshape(sigma2,[2,2]);
                
                    H=H+c_coupling(i1,i2,o)*Two_coupling(w(i2,1),w(i2,2),sigma1,sigma2);
                    o=o+1;
                
                end
            end
        end
    
        H_Matrix(i1,:,:)=H;
    end

    save(['Coe_',num2str(s),'.mat'],'c_self','c_coupling','c_time','H_Matrix')
end