function AIC=calculate_AIC(E,p)
    % calculating AIC coefficients
    [m,N,trials] = size(E);
    E_U=E;
    Sigma_prewave = zeros(m,m);
    for k = 1:size(E_U,3)
        Sigma(:,:,k) = E_U(:,:,k)*E_U(:,:,k)';
        Sigma_prewave = Sigma_prewave + Sigma(:,:,k);
    end
    Sigma_wave = Sigma_prewave./(N*k);
    Det_Sigma_wave = det(Sigma_wave);
    AIC = log(Det_Sigma_wave)+2*m*m*p/N;    
end
    