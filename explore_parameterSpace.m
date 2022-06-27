% load data in format channels by samples by trials
[m,N,k] = size(y);


%% estimate parameters
Uc=[0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
modelOrders=1:30;
lambda=Uc;
BIC_GLKF=zeros(size(modelOrders,2),size(Uc,2));
cons_GLKF=zeros(size(modelOrders,2),size(Uc,2));
siglev1_GLKF=zeros(size(modelOrders,2),size(Uc,2));
siglev2_GLKF=zeros(size(modelOrders,2),size(Uc,2));
stable_GLKF=zeros(size(modelOrders,2),size(Uc,2),340);
     

for i=1:size(modelOrders,2)
    p=modelOrders(i);
    for j=1:size(Uc,2)


    %% GLKF
    [ApAR, Awave, res]=GLKF_Eshwar(y,p, Uc(j));
    Awave=permute(Awave,[2,1,3]);

    res=permute(res,[2,3,1]);
    BIC_GLKF(i,j)=calculate_BIC(res,p);

    [cons_GLKF(i,j), ~] = consistency(y,ApAR);

    [siglev1_GLKF(i,j), siglev2_GLKF(i,j)]=arres_whiteness(res,p);

    stable_GLKF(i,j,:)=stability(Awave);

    %% RLS
%         [ApAR, Awave, res]=GLKF_RLS(y,p, lambda(j));
% 
%         res=permute(res,[2,3,1]);
%         BIC_RLS(i,j)=calculate_BIC(res,p);
% 
%         [cons_RLS(i,j), ~] = consistency(y,ApAR);
% 
%         [siglev1_RLS(i,j), siglev2_RLS(i,j)]=arres_whiteness(res,p);
% 
%         stable_RLS(i,j,:)=stability(Awave);

    end
end