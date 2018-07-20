function DTF_calculating(outputname,p)
    % Calculating DTF-parameters
    % Before using it, you should use 'GeneralLinearKalmanFilteringUni.m' to
    %   calculate all necessary parameters

    % GeneralLinearKalmanFilteringUni

    clearvars -except outputname p;  

    % Write what data you want to send in DTF { A - Abstractness_TimeSeries
    %                                           CmA - CommonAbstractness_TimeSeries
    %                                           CmC - CommonConcreteness_TimeSeries
    %                                           C - Concreteness_TimeSeries         }

    load(outputname)


    clear DTF FreqBrainMax f fmax fbin_per_Hz;

    Impone = eye( m*p , m*p );                                                 % make matrice of 1; dimension mp*mp
    Imone = eye( m , m );                                                      % make matrice of 1; dimension m*m
    Ikone = eye( k , k );                                                      % make matrice of 1; dimension k*k

    FreqBrainMax = 30;                                                         % MaxBrainFrequency - max freq of brain activity
    fmax = FreqBrainMax/freq;                                                  % maxf - max normalized frequency we need
    f = 0.0005:0.0005:fmax;
    fbin_per_Hz = size(f,2)/FreqBrainMax;                                      % bin_per_Hz - how many freq bins in 1 Hz
    for sn=1:N                                                                 % sn - sample number
        for lam=1:size(f,2)                                                    % lam - number of frequency bin
            clear prepA;                                                       %   lam is nonnormalized frequency; f is normalized
            prepA = 0;
            for q=1:p
                prepA = prepA + exp((-2i)*(pi)*(f(lam))*q) * ApAR(:,:,q,sn);
            end    
            AforH(:,:,lam,sn) = Imone - prepA;                                % AforH - matrix A to get matrix H
            H(:,:,lam,sn) = inv(AforH(:,:,lam,sn));                          % H - matrix we need for calclating DTF
            
            thetapre(:,:,lam,sn) = abs(H(:,:,lam,sn));
            theta(:,:,lam,sn) = (thetapre(:,:,lam,sn)).^2;
            
                                                                     % to channel
            for j=1:m                                                      % from channel
                DTF(:,j,lam, sn) = (theta(:,j,lam, sn)/(H(:,j,lam, sn)'*H(:,j,lam, sn));
            end
            
        end
    end

    % Here you can cave DTF Parameters

    
    save(outputname,'DTF','FreqBrainMax','f','fmax','fbin_per_Hz','-append');
    

end