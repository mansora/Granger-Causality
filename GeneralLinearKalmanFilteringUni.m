function GeneralLinearKalmanFilteringUni(s,p, outputname)
    %% Code created by Dmitriy Sakharuk
    % please contact me at "" for any comments or bug reporting
    % Calculating Autoregression parameters using General Linear Kalman Filtering
    % this code is based on the explanations provided in the paper 
    % "A Time-Varying Connectivity Analysis from Distributed EEG Sources: A Simulation Study" 
    % by Eshwar GHumare et al.
    
    % Input
    % s: input data in 3 dimensional .mat file
    % data format should be channels by time by trials
    % p: model order
    % outputname: name of the file with which the output will be save in
    % the same directory as the code
    
    % Output file will include

    % clear all nonuniversal values
    clearvars -except s p outputname;  
    

    

    % Write what data you want to send in GLKF { A - Abstractness_TimeSeries
    %                                            CmA - CommonAbstractness_TimeSeries
    %                                            CmC - CommonConcreteness_TimeSeries
    %                                            C - Concreteness_TimeSeries         }

    % dataname = 'A';
    
%     num_trials=[656, 614, 543, 625, 646, 295, 636, 650, 637, 546, 456, 656, 612];
    
    
    y = (10^12)*s;
%     y(1,1:550,:) = (10^10)*s(1:1,1:550,1:10);
%     y(2,1:550,:) = (10^10)*s(1:1,1:550,1:10);
%     for kk=1:10
%        y(:,1,kk) = zeros(2,1);
%        y(:,2,kk) = ones(2,1);
%     end
%     for i=1:500
%        for kk=1:10
%            y(:,i+2,kk) = [ (0.95*y(1,i+1,kk)) - (0.70*y(1,i,kk)) ;
%                            (0.50*y(1,i+1,kk)) - (0.90*y(1,i,kk)) ];
%        end            
%     end
    
%     for kk=1:10
%        y(:,1,kk) = zeros(2,1);
%        y(:,2,kk) = ones(2,1);
%     end
%     for i=1:500
%        for kk=1:10
%            y(:,i+2,kk) = [ (0.95*sqrt(2)*y(1,i+1,kk)) - (0.9025*y(1,i,kk)) ;
%                            (-0.50*y(2,i+1,kk)) + (0.5*y(1,i+1,kk)) ];
%        end            
%     end

    m = size(y,1);                                                             % m - number of channels
    % p = 2;                                                                     % p - GC model order
    N = size(y,2);                                                             % N - number of samples
    k = size(y,3);                                                             % k - number of trials
    time = 1.7;                                                                % time - timelength of data (in seconds)
    freq = N/time;                                                             % sampling frequency
%     dataname_for_save = dataname;                                              % just for saving dataname in file for future
    Uc = 0.02;                                                                 % The Update coefcient UC (0 < UC < 1) 
                                                                               %   controls the adaptation speed 
                                                                               %   of time-varying MVAR parameters A?p(n)
    
    Impone = eye( m*p , m*p );                                                 % make matrice of 1; dimension mp*mp
    Imone = eye( m , m );                                                      % make matrice of 1; dimension m*m
    Ikone = eye( k , k );                                                      % make matrice of 1; dimension k*k
   
                                                                                  %   all their values =0 or =1 where necessary
    Awave= zeros( m*p , m, N);                                                 % matrice of Kalman-Filtering MVAR parameters, =0;
                                                                               %   make start empty matrice Awave
    Hp = zeros( k , m*p, N);                                                   % make a matrice of Hp(1) coeff of 0.
                                                                               %   we won't use this matrice,
                                                                           %   but need smth in Hp(:,:,1);
                                                                           %   make start empty matrice Hp
    P = repmat(eye( m*p , m*p), [1 1 N]);                                           % a-posteriori error covariance matrix;    
                                                                           %   make start empty matrice P
    W= repmat(eye( m , m), [1 1 N]);                                               % the measurement error covariance;
                                                                           %   make start empty matrice W
    E= zeros( k , m, N);                                             % E - measurement error;
                                                                           %   make start empty matrice ErrorMeas
    X = zeros( k , k, N);                                             % the residual covariance;
                                                                           %   make start empty matrice X
    KG= zeros( m*p , k, N);                                          % the Kalman gain;
                                                                           %   make start empty matrice KG
    V= zeros( m*p , m*p, N);                                         % the state error covariance;
                                                                           %   make start empty matrice V
    %D(:,:,n) = ones( m , m );
    
    
   
 
    % place for D :D (because Debil)

    for n = (p+1):N
        
        for i = 1:p
            clear O;
            O = (squeeze(y(:,n-i,:)))';                                          % O - matrix of y for one sample n for all trials k 
            Hp(:,((i-1)*m+1):(i*m),n) = O;                                     % Hp - matrix of all previous values of y for order p
        end

%         for i = 1:p
%             clear O;
%             for L=1:k
%                 O(L,:) = y(:,n-i,L)';                                          % O - matrix of y for one sample n for all trials k
%             end    
%             Hp(:,((i-1)*m+1):(i*m),n) = O;                                     % Hp - matrix of all previous values of y for order p
%         end
%         for i=1:
%             yE(:,i) = y(:,n,i);                                                % yE - special matrix of y for calculating E(n)
%         end 
        
        yE = squeeze(y(:,n,:));                                                % yE - special matrix of y for calculating E(n)  
        E(:,:,n) = yE' - Hp(:,:,n)*Awave(:,:,n-1);                             % E - measurement error matrix
        W(:,:,n) = (1-Uc)*W(:,:,n-1) + Uc*(E(:,:,n)'*E(:,:,n))/(k-1) ;         % W - the measurement error covariance
        X(:,:,n) = inv(Hp(:,:,n)*P(:,:,n-1)*Hp(:,:,n)'+trace(W(:,:,n))*Ikone); % Calculate the residual covariance
        KG(:,:,n) = P(:,:,n-1)*Hp(:,:,n)'*X(:,:,n);                            %Calculate the Kalman gain
        Awave(:,:,n) = Awave(:,:,n-1)+KG(:,:,n)*E(:,:,n);                      %Update the MVAR estimates
        V(:,:,n) = ((Uc*trace((Impone-KG(:,:,n)*Hp(:,:,n))*P(:,:,n-1)))/(m*m*p))*Impone; %Calculate the state error covariance
        P(:,:,n) = (Impone-KG(:,:,n)*Hp(:,:,n))*P(:,:,n-1)+V(:,:,n);           %Update the a-posteriori error covariance
    end

    % Now we will transform parameters for AR-model from GLKF-parameters Awave

    ApAR(:,:,p,p) = zeros( m , m );

    for i=p+1:N
        for q=1:p
            ApAR(:,:,q,i) = Awave((m*(q-1)+1):(m*q),:,i)';                     % We take it from Awave, separating Awave by p
        end    
    end

    % Here you can save AutoRegression Parameters
    save(outputname,'ApAR','E','m','p','N','k','freq','time');

    
end
