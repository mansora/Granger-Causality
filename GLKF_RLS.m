function [ApAR, Awave, Z]=GLKF_RLS(y,p, lambda)
    %% recursive least squares algorithm for multitrial data, algorithm as specified in the paper of Milde et al (2010)
    % lambda=0.5;  % forgetting factor, needs to be value between 0 and 1
    % p - GC model order
    % m - number of channels
    % N - number of samples
    % k - number of trials

    [m,N,k] = size(y);
    Z=zeros(k,m,N);

    Awave = zeros( m , m*p, N);   
    C=zeros(m*p,m*p,N);
    W = zeros( k , m*p);   

    for n = (p+1):N
        
        for i = 1:p
            clear O;
            O = (squeeze(y(:,n-i,:)))';                                          % O - matrix of y for one sample n for all trials k 
            W(:,((i-1)*m+1):(i*m)) = O;                                     % Hp - matrix of all previous values of y for order p
        end
        C(:,:,n)=(1-lambda)*C(:,:,n-1)+W'*W;
        K=W/(C(:,:,n));
        Z(:,:,n)=squeeze(y(:,n,:))'-W*Awave(:,:,n-1)';
        Awave(:,:,n)=Awave(:,:,n-1)+Z(:,:,n)'*K;
        
    end

    % Now we will transform parameters for AR-model from GLKF-parameters Awave

    ApAR(:,:,p,p) = zeros( m , m );

    for i=p+1:N
        for q=1:p
            ApAR(:,:,q,i) = Awave(:,(m*(q-1)+1):(m*q),i)';                     % We take it from Awave, separating Awave by p
        end    
    end
    
    
end


