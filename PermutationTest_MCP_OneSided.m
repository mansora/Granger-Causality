function [SignificanceMap, SM]=PermutationTest_MCP_OneSided(Group1, Group2, Iterations, tthresh)
%% Permutation test with cluster based inference for multiple comparison for within and between subject analysis of time-frequency plots
    % TODO: make tthresh so that the user can select and if not default be
    % the following
    if nargin<4
        tthresh=2.2;
    end

    Group1=abs(Group1);
    Group2=abs(Group2);
    % TODO: add the unpaired ttest to possibilities
    [~,varr,fbin,tbin,n_pop]=size(Group1);
    % TODO: make sure size of group1 and group2 are the same
    Group=cat(5, Group1, Group2);
    
    % first compute actually observed test statistics
    STAT=@(a,b) ttest(a,b, 'Tail','right');
    Group1_temp=reshape(Group1,[],n_pop);
    Group2_temp=reshape(Group2,[],n_pop);
    [~,~,~,stats]=STAT(Group1_temp',Group2_temp');
    temp=stats.tstat;
    sigMap_Observed=reshape(temp,varr,varr,fbin,tbin);
    sigMap_Observed=double(sigMap_Observed>tthresh);
    
    
    MCP_obs=zeros(varr,varr);
    for i=1:varr
        for j=1:varr
            sigMap_Observed_temp=squeeze(sigMap_Observed(i,j,:,:));
            xtemp=bwconncomp(sigMap_Observed_temp);    
            clusterStats=zeros(1,xtemp.NumObjects);
            if xtemp.NumObjects~=0
                    for k=1:xtemp.NumObjects
                        clusterStats(k)=sum(sigMap_Observed_temp(xtemp.PixelIdxList{k}));
                    end
                    MCP_obs(i,j)=max(clusterStats); 
                else
                    MCP_obs(i,j)=max(max(sigMap_Observed_temp));
            end                       
        end
    end
    sigMap_Higher=zeros(varr,varr);
    sigMap_Lower=zeros(varr,varr);
    
    for It=1:Iterations
        x=randperm(n_pop*2);
        randPart_1=Group(:,:,:,:,x(1:n_pop));
        randPart_2=Group(:,:,:,:,x(n_pop+1:n_pop*2));
        
        randPart_1=reshape(randPart_1,[],n_pop);
        randPart_2=reshape(randPart_2,[],n_pop);
        [~,~,~,stats]=STAT(randPart_1',randPart_2');
        temp=stats.tstat;
        sigMap_temp=reshape(temp,varr,varr,fbin,tbin);
        
        sigMap_temp=double(sigMap_temp>tthresh);
        MCP_perm=zeros(varr,varr);
        for i=1:varr
            for j=1:varr
                sigMap_temp_temp=squeeze(sigMap_temp(i,j,:,:));
                xtemp=bwconncomp(sigMap_temp_temp); 
                clusterStats=zeros(1,xtemp.NumObjects);
                if xtemp.NumObjects~=0
                    for k=1:xtemp.NumObjects
                        clusterStats(k)=sum(sigMap_temp_temp(xtemp.PixelIdxList{k}));
                    end
                    MCP_perm(i,j)=max(clusterStats); 
                else
                    MCP_perm(i,j)=max(max(sigMap_temp_temp));
                end   
            end
        end
        
        sigMap_Higher=sigMap_Higher+double(MCP_perm>MCP_obs);
        sigMap_Lower=sigMap_Lower+double(MCP_perm<MCP_obs);
           
    end
    
    SignificanceMap=sigMap_Higher./sigMap_Lower;
    SM=zeros(varr,varr);
    SM(find(SignificanceMap<0.05))=1;
    
    figure,
    s1=0;
    clims=[0.1, 0.5];
    for i = 1 : varr
        for j = 1 : varr
            if i==j
               SM_tmp=zeros(1,1,fbin,tbin);      
            else
               SM_tmp=squeeze(SM(i,j,:,:));  
            end

            s1 = s1 + 1;
            h = subplot(varr,varr,s1);
            set(h,'FontSize',13);

            img = squeeze(SM_tmp);
            imagesc([0 30],[-0.2 1],img', clims) % On the frequency axis, 0.5 is equivalent with the Nyquist rate (Fs/2).
            colormap jet
            set(h,'YDir','normal')
            caxis([0 1])
            grid on

            if(i==varr && j==ceil(varr/2))
                xlabel('Normalized Frequency','Fontsize',14)
            end
            if(i==ceil(varr/2) && j==1)
                ylabel('Time (sample)','Fontsize',14)
            end
            title(['Ch' num2str(i) ' <--- Ch ' num2str(j)],'Fontsize',14)
        end
    end
    h = colorbar;
    set(h, 'Position', [.92 .11 .03 .8150],'FontSize',14)
end