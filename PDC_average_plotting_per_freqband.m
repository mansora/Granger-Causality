function PDC_average_plotting_per_freqband(outputname,i,j)%,subjects,dataname)
    % Plotting PDC for each frequence band for all time bins
    % Before using it, you should use 'GeneralLinearKalmanFilteringUni.m' to
    %   calculate all necessary parameters for AutoRegressive model
    % After that you should use 'PDC_calculating.m' to calculate all necessary
    %   parameters for PDC
    % And after that you should use 'PDC_average_per_freq.m' to calculate
    %   necessary average PDC parameters

    % GeneralLinearKalmanFilteringUni
    % PDC_calculating
    % PDC_average_per_freq

    clearvars -except outputname subjects dataname p i j;  

    % Write what data you want to send in PDC { A - Abstractness_TimeSeries
    %                                           CmA - CommonAbstractness_TimeSeries
    %                                           CmC - CommonConcreteness_TimeSeries
    %                                           C - Concreteness_TimeSeries         }

    load(outputname)

%     i = 2;   % to which channel
%     j = 1;   % from which channel
    plotsnums = 6;                                                             % number of plots = number of frequency bands
    x = linspace(-0.2,1.5,N);
    plotsnames = ['<PDC>f: channel y',int2str(j),'->y',num2str(i)];
%     figure('Name',plotsnames);
    for lam = 1:plotsnums                                
        clear z;
        if lam == 1
            for sn=1:N
                z(sn) = PDC_Average_Delta{sn}(i,j);
            end
        elseif lam == 2
            for sn=1:N
                z(sn) = PDC_Average_Theta{sn}(i,j);
            end   
        elseif lam == 3
            for sn=1:N
                z(sn) = PDC_Average_Alpha{sn}(i,j);
            end   
        elseif lam == 4
            for sn=1:N
                z(sn) = PDC_Average_Beta{sn}(i,j);
            end   
        elseif lam == 5
            for sn=1:N
                z(sn) = PDC_Average_Low_Gamma{sn}(i,j);
            end   
        elseif lam == 6
            for sn=1:N
                z(sn) = PDC_Average_High_Gamma{sn}(i,j);
            end 
        end    
        p = subplot(2,3,lam);
        plot(x,z);
        xlim([x(1) x(N)]);
        xlabel('Time, s');
        p.XGrid = 'on';
        ylim([0 1]);
        ylabel(['\langle\mid\pi_{',int2str(j),'\rightarrow',num2str(i),'}\mid^2\rangle_f']);
        if lam == 1
            title(['Delta band (0..4] Hz']);
        elseif lam == 2
            title(['Theta band (4..8] Hz']);
        elseif lam == 3
            title(['Alpha band (8..16] Hz']);
        elseif lam == 4
            title(['Beta band (16..32] Hz']);
        elseif lam == 5
            title(['Low Gamma band (32..48] Hz']);
        elseif lam == 6
            title(['High Gamma band (48..64] Hz']);
        end    
    end
%      figurename = ['Subject ',num2str(subjects),'. PDC Process ',dataname,'. Influence from ',num2str(j),' to ',num2str(i),' channels'];
% %     saveas(gcf,['Subject ',num2str(subjects),'. PDC Process ',dataname,'. Influence from ',num2str(j),' to ',num2str(i),' channels.png']);
%     print(figurename,'-dpng')
    
end

