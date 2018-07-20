function PDC_calculating(outputname,p)
    % Calculating PDC-parameters
    % Before using it, you should use 'GeneralLinearKalmanFilteringUni.m' to
    %   calculate all necessary parameters

    % GeneralLinearKalmanFilteringUni

    clearvars -except outputname p;  

    % Write what data you want to send in PDC { A - Abstractness_TimeSeries
    %                                           CmA - CommonAbstractness_TimeSeries
    %                                           CmC - CommonConcreteness_TimeSeries
    %                                           C - Concreteness_TimeSeries         }

    % dataname = 'C';
    load(outputname)
        
    clear PDC FreqBrainMax f fmax fbin_per_Hz;

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
            Adash{sn,lam}(:,:) = Imone - prepA;
            for i=1:m                                                          % to channel
                for j=1:m                                                      % from channel
                    prePDC{sn,lam}(i,j) = Adash{sn,lam}(i,j)/(sqrt(Adash{sn,lam}(:,j)'*Adash{sn,lam}(:,j)));
                    PDC{sn,lam}(i,j) = abs(prePDC{sn,lam}(i,j))^2;
                end
            end    
        end
    end

    % Here you can cave PDC Parameters

    
    save(outputname,'PDC','FreqBrainMax','f','fmax','fbin_per_Hz','-append');
    

end
