function [cons, reconData] = consistency(Y,ApAR)
    % Y=real data
    % ApAR estimated autoregressive parameters
    [m,N,k] = size(Y);
    p=size(ApAR,3);
%     %% compute residuals
%     for i=1:N
%         AR{i}=reshape(ApAR(:,:,:,i),[m],[]);
%     end
    
    add_noise=0;
    for trial=1:k
        sim_data = zeros(m,N);
        sim_data(:,1:p) = Y(:,1:p,trial);
    %     Err = zeros(m,N);

        for i=p+1:N
            for j=1:p
                sim_data(:,i)=sim_data(:,i)+ApAR(:,:,j,i)*sim_data(:,i-j);  
            end
            if add_noise
                sim_data(:,i) = awgn(sim_data(:,i), SNR, 'measured');
            end
        end
        X(:,:,trial)=sim_data;
    end
    
    
    s = k*(N-p);             % sample size (number of observations)


%     Y = X - E;                   % prediction
    reconData=X;
    X=X(:,:);
    Y=Y(:,:);
    
    Rr = (Y*Y')/(s-1);           % covariance estimate
    Rs = (X*X')/(s-1);           % covariance estimate

    cons = 1 - norm(Rs-Rr)/norm(Rr); % compare matrix norms


    

end



