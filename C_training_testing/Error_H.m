clear
clc

%% Generate basis, using BL represent basis_label
load('t_test.mat')

sigma_gather=zeros(4,2,2);
sigma_gather(1,:,:)=[1 0; 0 1];
sigma_gather(2,:,:)=[0 1; 1 0];
sigma_gather(3,:,:)=[0 -sqrt(-1); sqrt(-1) 0];
sigma_gather(4,:,:)=[1 0; 0 -1];

basis256=zeros(256,16,16);
BL=zeros(4,256);

o=1;
for i1=1:4
    
    x1=sigma_gather(i1,:,:);
    x1=reshape(x1,[2,2]);
    for i2=1:4
        
        x2=sigma_gather(i2,:,:);
        x2=reshape(x2,[2,2]);
        
        for i3=1:4
            
            x3=sigma_gather(i3,:,:);
            x3=reshape(x3,[2,2]);
            
            for i4=1:4
                
                x4=sigma_gather(i4,:,:);
                x4=reshape(x4,[2,2]);
                basis256(o,:,:)=kron(x1,kron(x2,kron(x3,x4)));
                BL(1,o)=i1;
                BL(2,o)=i2;
                BL(3,o)=i3;
                BL(4,o)=i4;
                o=o+1;
                
            end
        end
    end
end


for s=1:1
    
    eval(['load h_pred_',num2str(s),'.mat'])
    eval(['load Coe_',num2str(s),'.mat'])
    
    %% Change Hamiltonian to Schrodinger picture
    
    t=t_test;
    HH=h_pred;
    M=length(t);
    dt=t(2)-t(1);

    N=length(t);
    HS=zeros(N,256);
    HS(1,:)=HH(1,:);
    Ut=eye(16);
    for i=1:N-1
    
        HS_test=F_vector_matrix(HS(i,:));
    
        Ut=(eye(16)-sqrt(-1)*HS_test*dt)*Ut;
    
        HH_test=F_vector_matrix(HH(i,:));
        HS_test2=Ut*HH_test*Ut';
    
        HS(i+1,:)=F_matrix_vector(HS_test2);

    end

    HS=HS*5;

    T=zeros(256,256);
    HS_c=zeros(N,256);
    for i=1:256
        xb=basis256(i,:,:);
        xb=reshape(xb,[16,16]);
        f=F_matrix_vector(xb);
        T(:,i)=f;
    end
    
    %% Decouple
    for i=1:N
        % Tx=b
        b=HS(i,:);
        b=b';
        x=T\b;
        HS_c(i,:)=x;
    
    end
    
    %% self
    H_real(1,:)=c_self(1,1,1)+c_self(2,1,1)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_real(2,:)=c_self(1,1,2)+c_self(2,1,2)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_real(3,:)=c_self(1,1,3)+c_self(2,1,3)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    
	H_pred(1,:)=HS_c(:,65);
    H_pred(2,:)=HS_c(:,129);
    H_pred(3,:)=HS_c(:,193);
    
    %% coupling (first second)
    H_real(4,:)=c_coupling(1,1,1)+c_coupling(2,1,1)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(4,:)=HS_c(:,81);
    
	H_real(5,:)=c_coupling(1,1,2)+c_coupling(2,1,2)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(5,:)=HS_c(:,97);
    
	H_real(6,:)=c_coupling(1,1,3)+c_coupling(2,1,3)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(6,:)=HS_c(:,113);
    
	H_real(7,:)=c_coupling(1,1,4)+c_coupling(2,1,4)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(7,:)=HS_c(:,145);
    
	H_real(8,:)=c_coupling(1,1,5)+c_coupling(2,1,5)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(8,:)=HS_c(:,161);
    
	H_real(9,:)=c_coupling(1,1,6)+c_coupling(2,1,6)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(9,:)=HS_c(:,177);
    
	H_real(10,:)=c_coupling(1,1,7)+c_coupling(2,1,7)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(10,:)=HS_c(:,209);
    
	H_real(11,:)=c_coupling(1,1,8)+c_coupling(2,1,8)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(11,:)=HS_c(:,225);
    
	H_real(12,:)=c_coupling(1,1,9)+c_coupling(2,1,9)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(12,:)=HS_c(:,241);
    
	%% coupling (first fourth)
    H_real(13,:)=c_coupling(1,4,1)+c_coupling(2,4,1)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(13,:)=HS_c(:,66);
    
	H_real(14,:)=c_coupling(1,4,2)+c_coupling(2,4,2)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(14,:)=HS_c(:,67);
    
	H_real(15,:)=c_coupling(1,4,3)+c_coupling(2,4,3)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(15,:)=HS_c(:,68);
    
	H_real(16,:)=c_coupling(1,4,4)+c_coupling(2,4,4)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(16,:)=HS_c(:,130);
    
	H_real(17,:)=c_coupling(1,4,5)+c_coupling(2,4,5)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(17,:)=HS_c(:,131);
    
	H_real(18,:)=c_coupling(1,4,6)+c_coupling(2,4,6)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(18,:)=HS_c(:,132);
    
	H_real(19,:)=c_coupling(1,4,7)+c_coupling(2,4,7)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(19,:)=HS_c(:,194);
    
	H_real(20,:)=c_coupling(1,4,8)+c_coupling(2,4,8)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(20,:)=HS_c(:,195);
    
	H_real(21,:)=c_coupling(1,4,9)+c_coupling(2,4,9)*sin(c_time(2)*t(1:N)+2*pi*c_time(3));
    H_pred(21,:)=HS_c(:,196);
    
    %% coupling (three body)
    o=22;
    for i=1:256
        count_BL_larger1=0;
        for j=1:4
            if BL(j,i)>1
                count_BL_larger1=count_BL_larger1+1;
            end
        end
        if BL(1,i)~=1 & count_BL_larger1>2
            H_pred(o,:)=HS_c(:,i);
            H_real(o,:)=0;
            o=o+1;
        end
        if BL(1,i)~=1 & BL(2,i)==1 & BL(3,i)~=1 & BL(4,i)==1
            H_pred(o,:)=HS_c(:,i);
            H_real(o,:)=0;
            o=o+1;
        end
    end
    
    %H_pred=-H_pred;
    error_Fo_prime(s)=mean(mean(abs(H_real-H_pred).^2))/mean(mean(abs(H_real).^2));
    disp([s,mean(error_Fo_prime)])
    
%     %% plot
%     for i=1:21
%         figure()
%         plot(t,H_real(i,:),'linewidth',1.5);hold on
%         plot(t,H_pred(i,:),'--','linewidth',1.5);
%     end
%     figure()
%     for i=22:length(H_real(:,1))
%         plot(t,H_pred(i,:),'--','linewidth',1.5);hold on
%     end
%     
end

