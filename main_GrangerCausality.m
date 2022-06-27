clear all
close all
load('E:\Algorithm for Mansoureh\GrangerRLS\Data\TimeSeriesFinalAllExtreme.mat')
for subb=1:21
p=14;
clearvars -except subb AbstractnessPerSubject ConcretenessPerSubject
% subb=1;

%% load data

% load('F:\Raw Data Categorization and NaturalReading\preprocessedCategorizationTask\data150\results\TimeSeriesCreation\TimeSeriesCommonAERPreduced.mat')
% y=cat(3,AbstractnessPerSubject{subb},ConcretenessPerSubject{subb});

% load('F:\Raw Data Categorization and NaturalReading\preprocessedCategorizationTask\data150\results\TimeSeriesCreation\TimeSeriesCommon.mat')
% load('E:\Algorithm for Mansoureh\GrangerRLS\Data\TimeSeriesCommonExtreme.mat')
% load('E:\Algorithm for Mansoureh\GrangerRLS\Data\TimeSeriesRandomExtreme.mat')
% load('E:\Algorithm for Mansoureh\GrangerRLS\Data\TimeSeriesFinalAll.mat')


% % load('F:\Raw Data Categorization and NaturalReading\preprocessedCategorizationTask\data150\results\TimeSeriesCreation\TimeSeriesCommonMean.mat')
% y=cat(3,AbstractnessPerSubject{subb},ConcretenessPerSubject{subb});
y=AbstractnessPerSubject{subb};
y=y(1:8,:,:);
% % y2=cat(3,AbstractnessPerSubject{subb+1},ConcretenessPerSubject{subb+1});
% % y=cat(3,y1,y2);
y=y*(10^10);


% load('F:\Raw Data Categorization and NaturalReading\preprocessedCategorizationTask\data150\Data_Subject1.mat')
% load('F:\Raw Data Categorization and NaturalReading\preprocessedCategorizationTask\channel_BrainProducts_EasyCap128.mat')
% load('F:\Raw Data Categorization and NaturalReading\preprocessedCategorizationTask\data150\results\TimeSeriesCreation\AbConMarker.mat')
% y=data([49,53,65,71,87,92],:,:);
% MarkerSub1=AbConMarker(1:656);
% y=y(:,:,MarkerSub1~=0);

[m,N,k] = size(y);

%% estimate parameters

% lambda=0.2;
% [ApAR, Awave, res]=GLKF_RLS(y,p, lambda);

% E=0.7;
% Uc=0.02;
% [ApAR, Awave]=GLKF_Milde_DS(y,p, E, Uc);
 
Uc=0.001;
[ApAR, Awave, res]=GLKF_Eshwar(y,p, Uc);
Awave=permute(Awave,[2,1,3]);


%% model order validation using Bayesian Information Criterion


% res=permute(res,[2,3,1]);
% BIC=calculate_BIC(res,p);

%% whiteness test
% [siglev_Mode1, siglev_Mode2]=arres_whiteness(res,p);
% siglev_Mode1/20


%% Consistency test

% [cons, X] = consistency(y,ApAR);


%% should be stability test, for now we plot 1 random trial to see stability

% plot_one_simulated_trial(100,y,p,ApAR)


%% PDC and DTF

p_opt=p;
Fs=200;
Fmax=Fs/6;
Nf=100;  % number of frequency bins
CH=m;



Down_Sample=1;
Awave = Awave(:,:,1:Down_Sample:end); % Down sampling in the AR coefficients
T = size(Awave,3);

[PDC_TV, DTF_TV] = PDC_DTF_matrix(Awave,p_opt,Fs,Fmax,Nf); % for now only works with GLKF_RLS


%% Surrogate data method (shuffling over samples)
N_Surr = 1; % Number of surrogates
PDC_Surr = zeros([size(PDC_TV) N_Surr]); % CH x CH x N_freq x N_segments x N_surr
DTF_Surr = zeros([size(DTF_TV) N_Surr]); % CH x CH x N_freq x N_segments x N_surr

for surr = 1 : N_Surr
    
    clear A_tmp A_tmp_reshape y_surr
%     for trial=1:k
%         for j = 1 : CH
%             y_surr(j,:,trial) = y(j,randperm(size(y,2)),trial); % Randomize all columns of 'y'
%         end
%     end
    
    % creating surrogate data using phase randomization method from Tim
    % Mullen (SIFT toolbox)
    for trial=1:k
        y_surr(:,:,trial) = ...
        ifft(abs(fft(y(:,:,trial))) ...
        .* exp(1i*angle(fft(rand(m,N)))), ...
        'symmetric');
    end

    [~, A_tmp, ~]=GLKF_Eshwar(y_surr,p, Uc);
    A_tmp=permute(A_tmp,[2,1,3]);
    A_tmp2 = A_tmp(:,:,1:Down_Sample:end);

    [PDC_Surr(:,:,:,:,surr), DTF_Surr(:,:,:,:,surr)] = PDC_DTF_matrix(A_tmp2,p_opt,Fs,Fmax,Nf);
    
    disp(['The surrogate number ' num2str(surr) ' was finished.'])

end

PDC_Surr2 = max(PDC_Surr,[],5); % Maximum values of the PDC measure at each time and each frequency over surrogates for all channels are extracted.
DTF_Surr2 = max(DTF_Surr,[],5); % Maximum values of the DTF measure at each time and each frequency over surrogates for all channels are extracted.

%% Plot - Adaptive PDC and DTF (time-varying)
figure, % ---> DTF plot 
s1 = 0;
alpha_level = .01; % Alpha level (1 - confidence interval) ---> for hypothesis testing
L=N;
clims=[0.1, 0.5];
for i = 1 : CH
    for j = 1 : CH
        if i==j
           DTF_tmp=zeros(1,1,Nf,T);
        else
            DTF_tmp = abs(DTF_TV(i,j,:,:));
            DTF_Surr_tmp = abs(DTF_Surr2(i,j,:,:));
            mask = DTF_tmp>(1-alpha_level/2)*DTF_Surr_tmp;
            DTF_tmp=DTF_tmp.*mask;
        end
        
        s1 = s1 + 1;
        h = subplot(CH,CH,s1);
        set(h,'FontSize',13);
        
        img = squeeze(DTF_tmp);
        imagesc([0 30],[-0.2 1],img',clims) % On the frequency axis, 0.5 is equivalent with the Nyquist rate (Fs/2).
        colormap jet
        set(h,'YDir','normal')
        caxis([0 1])
        grid on
        
        if(i==CH && j==ceil(CH/2))
            xlabel('Normalized Frequency','Fontsize',14)
        end
        if(i==ceil(CH/2) && j==1)
            ylabel('Time (sample)','Fontsize',14)
        end
        title(['Ch' num2str(i) ' <--- Ch ' num2str(j)],'Fontsize',14)
    end
end
h = colorbar;
set(h, 'Position', [.92 .11 .03 .8150],'FontSize',14)

filenamePlotDTF1=['Concrete_DTF_Subject',num2str(subb),'order',num2str(p),  '.png'];
filenamePlotDTF2=['Concrete_DTF_Subject',num2str(subb),'order',num2str(p),'.fig'];
saveas(gcf, filenamePlotDTF1)
saveas(gcf, filenamePlotDTF2)



figure, % ---> PDC plot 
s1 = 0;
clear mask
for i = 1 : CH
    for j = 1 : CH
        if i==j
           PDC_tmp=zeros(1,1,Nf,T);
        else
            
            PDC_tmp = abs(PDC_TV(i,j,:,:));
            PDC_Surr_tmp = abs(PDC_Surr2(i,j,:,:));
            mask = PDC_tmp>(1-alpha_level/2)*PDC_Surr_tmp;
            PDC_tmp=PDC_tmp.*mask;
        end
        
        s1 = s1 + 1;
        h = subplot(CH,CH,s1);
        set(h,'FontSize',13);
        
        img = squeeze(PDC_tmp);
        imagesc([0 30],[-0.2 1],img', clims) % On the frequency axis, 0.5 is equivalent with the Nyquist rate (Fs/2).
        colormap jet
        set(h,'YDir','normal')
        caxis([0 1])
        grid on
        
        if(i==CH && j==ceil(CH/2))
            xlabel('Normalized Frequency','Fontsize',14)
        end
        if(i==ceil(CH/2) && j==1)
            ylabel('Time (sample)','Fontsize',14)
        end
        title(['Ch' num2str(i) ' <--- Ch ' num2str(j)],'Fontsize',14)
    end
end
h = colorbar;
set(h, 'Position', [.92 .11 .03 .8150],'FontSize',14)

filenamePlotPDC1=['Concrete_PDC_Subject',num2str(subb),'order',num2str(p),  '.png'];
filenamePlotPDC2=['Concrete_PDC_Subject',num2str(subb),'order',num2str(p),'.fig'];
saveas(gcf, filenamePlotPDC1)
saveas(gcf, filenamePlotPDC2)

end





