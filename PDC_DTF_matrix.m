function [PDC, DTF] = PDC_DTF_matrix(A,p_opt,Fs,Fmax,Nf)
% This function computes time-varying PDC and DTF measures based on the
% frozen-time MVAR models.
%
% Reference paper: A. Omidvarnia, M. Mesbah, J. M. O'Toole et al.,
% �Analysis of the time-varying cortical neural connectivity in the newborn EEG: A time-frequency approach,�
% in Systems, Signal Processing and their Applications (WOSSPA), 2011 7th International Workshop on, 2011, pp. 179-182
%
% The paper is also available at:
% http://qspace.qu.edu.qa/handle/10576/10737?show=full
%
% **** Input arguments:
% A (CH x CH*p): MVAR matrices associated with delays from r=1 to r=p. 
% 'CH' is the number of channels and A = [A1 A2 ... Ap].
%
% p_opt: Optimum MVAR model order which can be determined by AIC or BIC measures.
%
% Fs: Sampling frequency
%
% Nf: Number of frequency bins
%
% Fmax: maximum frequency limit (must be less than Fs/2)
%
% **** Output arguments:
% PDC (CH x CH x Nf x T) : a 4-D array containing the PDC values of all channel combinations at all time-frequency bins.
%
% DTF (CH x CH x Nf x T) : a 4-D array containing the DTF values of all channel combinations at all time-frequency bins.
%
%
% Written by: Amir Omidvarnia
%
% See also: BioSig toolbox (available at: http://biosig.sourceforge.net/).

[CH,L,T] = size(A);

PDC   = zeros(CH,CH,Nf,T); % Partial Directed Coherence
DTF   = zeros(CH,CH,Nf,T); % Directed Transfer Function

f = (0:Nf-1)*(Fmax/Nf); % Frequency span
z = 1i*2*pi/Fs;

for t = 1 : T % Number of time points
    
    A2 = [eye(CH) -squeeze(A(:,:,t))]; % A2(f) = [I -A1(t) -A2(t) ... -Ap(t)]

    %% **** Partial Directed Coherence (PDC)
    for n = 1 : Nf % Number of frequency points

        % **** Matrix A(f) --> Summation of the AR coefficients matrices in the frequency domain
        Af = zeros(CH);
        for k = 1 : p_opt+1
            Af = Af + A2(:,k*CH+(1-CH:0))*exp(z*(k-1)*f(n));
        end

        %% **** Matrix H(f) --> Transfer Function matrix (necessary for computing DTF)
        H(:,:,n) = inv(Af);

        %% **** Partial Directed Coherence function
        for ch = 1 : CH
            tmp = squeeze(Af(:,ch)); % ii'th column of A(f) --> a_j
            denom_PDC(ch) = sqrt(tmp'*tmp); % Hermitian(a_j) * a_j, j = 1:CH --> denominator of the PDC function
        end;

        PDC(:,:,n,t)  = abs(Af)./denom_PDC(ones(1,CH),:); % Aij/(sqrt(aj'*aj))

    end

    %% **** Directed Transfer Function (DTF)
    for chi = 1 : CH

        denum_DTF  = sum(abs(H(chi,:,:)).^2,2); % (Hi1^2+Hi2^2+...)

        for chj = 1 : CH
            DTF(chi,chj,:,t)   = abs(H(chi,chj,:))./sqrt(denum_DTF); % Directed Transfer Function: Hij / sqrt(Hi1^2+Hi2^2+...)
        end

    end
    %% ****
    
end


