% close all
L=1000;
y=zeros(3,L,100);
 for kk=1:10
   y(:,1,kk) = zeros(3,1);
   y(:,2,kk) = ones(3,1);
 end

 
     
 for kk=1:100
     w=rand(L,3);
    mean_w=mean(w,1);
     w(:,1)=w(:,1)-mean_w(1);
     w(:,2)=w(:,2)-mean_w(2);
     w(:,3)=w(:,3)-mean_w(3);
    for n=3:L
       if n<L/2
           c12(n)=0.5*(n/L);
       else
           c12(n)=0.5*(L-n)/(L/2);
       end

       if n<0.7*L
           c23(n)=0.4;
       else
           c23(n)=0;
       end

       y(:,n,kk) = [ 0.5*y(1,n-1,kk)  - 0.7*y(1,n-2,kk) + c12(n)*y(2,n-1,kk) + w(n,1);
                     0.7*y(2,n-1,kk)  - 0.5*y(2,n-2,kk) + 0.2*y(1,n-1,kk) + c23(n)*y(3,n-1,kk) + w(n,2);
                     0.8*y(3,n-1,kk)  + w(n,3)];                 
    end 
 end

    figure,
    for i=1:100
        subplot(311), plot(y(1,:,i))
        hold on
        subplot(312), plot(y(2,:,i))
        hold on
        subplot(313), plot(y(3,:,i))
        hold on
    end
        
        



% for kk=1:10
%        y(:,1,kk) = zeros(2,1);
%        y(:,2,kk) = ones(2,1);
%     end
%     for i=1:500
%        for kk=1:10
%            y(:,i+2,kk) = [ (0.95*y(1,i+1,kk)) - (0.70*y(1,i,kk)) ;
%                            (0.50*y(1,i+1,kk)) - (0.90*y(1,i,kk)) ];
%        end            
%     end
%     
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
% 
% %     y=zeros(3,502,10);
% %      for kk=1:10
% %        y(:,1,kk) = zeros(3,1);
% %        y(:,2,kk) = ones(3,1);
% %      end
% %      w=rand(502,3);
% %      mean_w=mean(w,1);
% %      w(:,1)=w(:,1)-mean_w(1);
% %      w(:,2)=w(:,1)-mean_w(2);
% %      w(:,3)=w(:,1)-mean_w(3);
% %     
% %     for i=2:500
% %        for kk=1:10
% %            y(:,i+1,kk) = [ (0.5*y(1,i,kk))  + (0.3*y(2,i,kk)) + (0.4*y(3,i,kk)) + w(i,1);
% %                            (-0.5*y(1,i,kk)) + (0.3*y(2,i,kk)) + (1*y(3,i,kk))   + w(i,2);
% %                            (-0.3*y(1,i,kk)) - (0.2*y(2,i,kk)) + (1*y(3,i,kk))   + w(i,3)];
% %        end            
% %     end

% clear all
% close all
p=15;
subb=10;

% load('F:\Raw Data Categorization and NaturalReading\preprocessedCategorizationTask\data150\results\TimeSeriesCreation\TimeSeriesCommonAERPreduced.mat')
% y=cat(3,AbstractnessPerSubject{subb},ConcretenessPerSubject{subb});

% load('F:\Raw Data Categorization and NaturalReading\preprocessedCategorizationTask\data150\results\TimeSeriesCreation\TimeSeriesCommon.mat')
% y=cat(3,AbstractnessPerSubject{subb},ConcretenessPerSubject{subb});
% y2=cat(3,AbstractnessPerSubject{subb+1},ConcretenessPerSubject{subb+1});
% y=cat(3,y1,y2);

% y=y*(10^10);
y(1,:,:)=-y(1,:,:);
[m,N,k] = size(y);

%% estimate parameters

% lambda=0.01;
% [ApAR, Awave]=GLKF_RLS(y,p, lambda);

% E=0.7;
% Uc=0.02;
% [ApAR, Awave]=GLKF_Milde_DS(y,p, E, Uc);
 
Uc=0.01;
[ApAR, Awave, res]=GLKF_Eshwar(y,p, Uc);
Awave=permute(Awave,[2,1,3]);

%% model order validation using Bayesian Information Criterion

% for i=1:N
%     AR{i}=reshape(ApAR(:,:,:,i),[m],[]);
% end
% res=est_mvarResiduals(y,ApAR);
% BIC=calculate_BIC(res,p)

%% Consistency test

% [cons, X] = consistency(y,ApAR)


%% should be stability test, for now we plot 1 random trial to see stability
% plot_one_simulated_trial(300,y,p,ApAR)


%% PDC and DTF

p_opt=p;
Fs=200;
Fmax=Fs/2;
Nf=30;  % number of frequency bins
CH=m;
Down_Sample=10;
Awave = Awave(:,:,1:Down_Sample:end); % Down sampling in the AR coefficients
T = size(Awave,3);

[PDC_TV, DTF_TV] = PDC_DTF_matrix(Awave,p_opt,Fs,Fmax,Nf); % for now only works with GLKF_RLS


%% Surrogate data method (shuffling over samples)
N_Surr = 10; % Number of surrogates
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

    PDC_Surr2 = max(PDC_Surr,[],5);


%% plot PDC and DTF

% % Plot - Adaptive PDC and DTF (time-varying)
% figure, % ---> DTF plot 
% s1 = 0;
% alpha_level = .01; % Alpha level (1 - confidence interval) ---> for hypothesis testing
% L= N; %size(y2,2);
% for i = 1 : CH
%     for j = 1 : CH
%         s1 = s1 + 1;
%         h = subplot(CH,CH,s1);
%         set(h,'FontSize',13);
%         
%         DTF_tmp = abs(DTF_TV(i,j,:,:));
% %         DTF_Surr_tmp = abs(DTF_Surr2(i,j,:,:));
% %         mask = DTF_tmp>(1-alpha_level/2)*DTF_Surr_tmp;
%         
%         img = squeeze(DTF_tmp);
%         imagesc([0 .17],[0 L],img') % On the frequency axis, 0.5 is equivalent with the Nyquist rate (Fs/2).
%         colormap jet
%         set(h,'YDir','normal')
%         caxis([0 1])
%         grid on
%         
%         if(i==CH && j==ceil(CH/2))
%             xlabel('Normalized Frequency','Fontsize',14)
%         end
%         if(i==ceil(CH/2) && j==1)
%             ylabel('Time (sample)','Fontsize',14)
%         end
%         title(['Ch' num2str(i) ' <--- Ch ' num2str(j)],'Fontsize',14)
%     end
% end
% h = colorbar;
% set(h, 'Position', [.92 .11 .03 .8150],'FontSize',14)
alpha_level=0.05;
figure, % ---> PDC plot 
s1 = 0;
L=N;
clear mask
for i = 1 : CH
    for j = 1 : CH
        if i==j
               PDC_tmp=zeros(1,1,30,100);      
            else
               PDC_tmp = abs(PDC_TV(i,j,:,:));
        end
            
        s1 = s1 + 1;
        h = subplot(CH,CH,s1);
        set(h,'FontSize',13);
        
        
        PDC_Surr_tmp = abs(PDC_Surr2(i,j,:,:));
        mask = PDC_tmp>(1-alpha_level/2)*PDC_Surr_tmp;
        
        img = squeeze(PDC_tmp);
        imagesc([0 1],[0 L],img') % On the frequency axis, 0.5 is equivalent with the Nyquist rate (Fs/2).
        colormap jet
        set(h,'YDir','normal')
        caxis([0 1])
        grid on
        im=gca;
        im.FontSize=16;
        im.FontWeight='bold';
        im.LineWidth=2;
        
        if(i==CH && j==ceil(CH/2))
            xlabel('Normalized Frequency','Fontsize',16, 'Fontweight', 'bold')
        end
        if(i==ceil(CH/2) && j==1)
            ylabel('Time (sample)','Fontsize',16, 'Fontweight', 'bold')
        end
        title(['variable' num2str(i) ' <- variable ' num2str(j)],'Fontsize',16, 'Fontweight', 'bold')
    end
end
h = colorbar;
set(h, 'Position', [.92 .11 .03 .8150],'FontSize',16, 'Fontweight', 'bold',  'LineWidth', 2)


    
    
 
