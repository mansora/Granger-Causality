function BIC=calculate_BIC(E,p)
    % Calculating BIC-coefficients   
    [m,N,trials] = size(E);
    E_U = E;  %permute(E,[2,3,1]);
    Sigma_prewave = zeros(m,m);
    for k = 1:size(E_U,3)
        Sigma(:,:,k) = E_U(:,:,k)*E_U(:,:,k)';
        Sigma_prewave = Sigma_prewave + Sigma(:,:,k);
    end
    Sigma_wave = Sigma_prewave./(N*k);
    Det_Sigma_wave = det(Sigma_wave);
    BIC = log(Det_Sigma_wave)+log(N)*m*m*p/N;
    
end
    