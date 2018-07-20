%                                                  A - Abstractness_TimeSeries
    %                                            CmA - CommonAbstractness_TimeSeries
    %                                            CmC - CommonConcreteness_TimeSeries
    %                                            C - Concreteness_TimeSeries         }

%% Input your model order p
CurrDir=pwd;
p = 5;
% data format should be channels by time by trials

%%   
% data_is_downloaded = exist('Abstractness_TimeSeries');
% if data_is_downloaded == 0
%     load('FinalTimeSeries_Data_withDescription.mat');
% end
figure;
num_trials=10;

list={'Abstract', 'CommAbstract', 'CommConcrete', 'Concrete'};
for nameprocess=1:4
    dataname=list{nameprocess};
    sTotal=[];
    outputDir = sprintf('%s\\First10TrialsPooled\\%s\\',CurrDir, dataname);
    % Check if the folder exists, and if not, make it...
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

        
           
        if strcmp(dataname, 'Abstract')
%             k=strfind(Abstractness_Description,['subject', num2str(subjects),'/']);
%             kk=find(~cellfun(@isempty,k));
%             s = permute(Abstractness_TimeSeries,[3 2 1]);
%             s=s(:,:,kk); % take all trials in subject
%             s=s(:,:,kk(1:num_trials)); % take only first #num_trials of trials in each subject
%             sTotal=cat(3,sTotal,s);
            outputname=[outputDir, 'Abstractness_Results'];
            m=9;
        elseif strcmp(dataname, 'CommAbstract')
%             k=strfind(CommonAbstractness_Description,['subject', num2str(subjects),'/']);
%             kk=find(~cellfun(@isempty,k));
%             s = permute(CommonAbstractness_TimeSeries,[3 2 1]); 
%             s=s(:,:,kk); % take all trials in subject
%             s=s(:,:,kk(1:num_trials)); % take only first #num_trials of trials in each subject
%             sTotal=cat(3,sTotal,s);
            outputname=[outputDir, 'CommonAbstractness_Results'];
            m=4;
        elseif strcmp(dataname, 'CommConcrete')
%             k=strfind(CommonConcreteness_Description,['subject', num2str(subjects),'/']);
%             kk=find(~cellfun(@isempty,k));
%             s = permute(CommonConcreteness_TimeSeries,[3 2 1]);
%             s=s(:,:,kk); % take all trials in subject
%             s=s(:,:,kk(1:num_trials)); % take only first #num_trials of trials in each subject
%             sTotal=cat(3,sTotal,s);
            outputname=[outputDir, 'CommonConcreteness_Results'];
            m=4;
        elseif strcmp(dataname, 'Concrete')
%             k=strfind(Concreteness_Description,['subject', num2str(subjects),'/']);
%             kk=find(~cellfun(@isempty,k));
%             s = permute(Concreteness_TimeSeries,[3 2 1]);
%             s=s(:,:,kk); % take all trials in subject
%             s=s(:,:,kk(1:num_trials)); % take only first #num_trials of trials in each subject
%             sTotal=cat(3,sTotal,s);
            outputname=[outputDir, 'Concreteness_Results'];
            m=8;
        end
        
        plottingTF(outputname)
            
    
    
end



