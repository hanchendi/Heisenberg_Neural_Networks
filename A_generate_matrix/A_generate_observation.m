clear
clc

sigma0=[1 0; 0 1];
sigma1=[0 1; 1 0];
sigma2=[0 -sqrt(-1); sqrt(-1) 0];
sigma3=[1 0; 0 -1];

basis4=zeros(3,16,16);

basis4(1,:,:)=kron(sigma1,eye(8));
basis4(2,:,:)=kron(sigma2,eye(8));
basis4(3,:,:)=kron(sigma3,eye(8));

save operater_observe.mat basis4