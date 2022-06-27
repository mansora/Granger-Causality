function [ApAR, Awave, E]=GLKF_Eshwar(y,p, Uc)
    %% GENERAL LINEAR KALMAN FILTER for multitrial data, algorithm as specified in the paper of Ghumare et al (2017)
    % p  - GC model order
    % Uc - The Update coefcient UC (0 < UC < 1) controls the adaptation
    % speed of time-varying parameters
    % m  - number of channels
    % N  - number of samples
    % k  - number of trials
    [m,N,k] = size(y);
  
    Impone = eye( m*p , m*p );                                                 % make matrice of 1; dimension mp*mp
    Ikone = eye( k , k );                                                      % make matrice of 1; dimension k*k
   
                                                                                  %   all their values =0 or =1 where necessary
    Awave = zeros( m*p , m, N);                                                 % matrice of Kalman-Filtering MVAR parameters, =0;
                                                                               %   make start empty matrice Awave
    Hp = zeros( k , m*p, N);                                                   % make a matrice of Hp(1) coeff of 0.
                                                                               %   we won't use this matrice,
                                                                           %   but need smth in Hp(:,:,1);
                                                                           %   make start empty matrice Hp
    P = repmat(eye( m*p , m*p), [1 1 N]);                                           % a-posteriori error covariance matrix;    
                                                                           %   make start empty matrice P
    W = repmat(eye( m , m), [1 1 N]);                                               % the measurement error covariance;
                                                                           %   make start empty matrice W
    E = zeros( k , m, N);                                             % E - measurement error;
                                                                           %   make start empty matrice ErrorMeas
    X = zeros( k , k, N);                                             % the residual covariance;
                                                                           %   make start empty matrice X
    KG = zeros( m*p , k, N);                                          % the Kalman gain;
                                                                           %   make start empty matrice KG
    V = zeros( m*p , m*p, N);                                         % the state error covariance;
                                                                           %   make start empty matrice V
    
   
    for n = (p+1):N
        
        for i = 1:p
            clear O;
            O = (squeeze(y(:,n-i,:)))';                                          % O - matrix of y for one sample n for all trials k 
            Hp(:,((i-1)*m+1):(i*m),n) = O;                                     % Hp - matrix of all previous values of y for order p
        end


        
        yE = squeeze(permute(y(:,n,:),[1,3,2]));                               % yE - special matrix of y for calculating E(n)  
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
            ApAR(:,:,q,i) = Awave((m*(q-1)+1):(m*q),:,i);                     % We take it from Awave, separating Awave by p
        end    
    end
    
end