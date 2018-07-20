%                                                  A - Abstractness_TimeSeries
    %                                            CmA - CommonAbstractness_TimeSeries
    %                                            CmC - CommonConcreteness_TimeSeries
    %                                            C - Concreteness_TimeSeries         }

%% Input your model order p
CurrDir=pwd;
p = 8;
% data format should channels by time by trials

%%   
data_is_downloaded = exist('Abstractness_TimeSeries');
if data_is_downloaded == 0
    load('FinalTimeSeries_Data_withDescription.mat');
end

num_trials=20;
list={'Abstract', 'CommAbstract', 'CommConcrete', 'Concrete'};
for subjects=1:13  % 1:13
    outputDir = sprintf('%s\\AllTrialsPerSubject\\Subject%s\\',CurrDir, num2str(subjects));
    % Check if the folder exists, and if not, make it...
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    

    for i=1:4
        dataname=list{i};
        if strcmp(dataname, 'Abstract')
            k=strfind(Abstractness_Description,['subject', num2str(subjects),'/']);
            kk=find(~cellfun(@isempty,k));
            s = permute(Abstractness_TimeSeries,[3 2 1]);
            s=s(:,:,kk); % take all trials in subject
%             s=s(:,:,kk(1:num_trials)); % take only first #num_trials of trials in each subject
            outputname=[outputDir, 'Abstractness_Results'];
        elseif strcmp(dataname, 'CommAbstract')
            k=strfind(CommonAbstractness_Description,['subject', num2str(subjects),'/']);
            kk=find(~cellfun(@isempty,k));
            s = permute(CommonAbstractness_TimeSeries,[3 2 1]); 
            s=s(:,:,kk); % take all trials in subject
%             s=s(:,:,kk(1:num_trials)); % take only first #num_trials of trials in each subject
             outputname=[outputDir, 'CommonAbstractness_Results'];
        elseif strcmp(dataname, 'CommConcrete')
            k=strfind(CommonConcreteness_Description,['subject', num2str(subjects),'/']);
            kk=find(~cellfun(@isempty,k));
            s = permute(CommonConcreteness_TimeSeries,[3 2 1]);
            s=s(:,:,kk); % take all trials in subject
%             s=s(:,:,kk(1:num_trials)); % take only first #num_trials of trials in each subject
            outputname=[outputDir, 'CommonConcreteness_Results'];
        elseif strcmp(dataname,  'Concrete')
            k=strfind(Concreteness_Description,['subject', num2str(subjects),'/']);
            kk=find(~cellfun(@isempty,k));
            s = permute(Concreteness_TimeSeries,[3 2 1]);
            s=s(:,:,kk); % take all trials in subject
%             s=s(:,:,kk(1:num_trials)); % take only first #num_trials of trials in each subject
            outputname=[outputDir, 'Concreteness_Results'];
        end


        fprintf('General Linear Kalman Filtering for %s \n', outputname);
        GeneralLinearKalmanFilteringUni(s, p, outputname)
        fprintf('PDC calculating for %s \n', outputname);
        PDC_calculating(outputname,p)
        PDC_average_per_freq(outputname)
%         PDC_average_plotting_per_freqband(outputname)
        
        fprintf('DTF calculating for %s \n', outputname);
        DTF_calculating(outputname,p)
        DTF_average_per_freq(outputname)
%         DTF_average_plotting_per_freqband(outputname)
    end

end




