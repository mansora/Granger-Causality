function plot_granger(SM)
[~,varr, fbin,tbin_n_pop]=size(SM);
%% just plot significance
figure,
    s1=0;
%     clims=[0.1, 0.5];
    for i = 1 : varr
        for j = 1 : varr
            if i==j
               SM_tmp=zeros(1,1,fbin,tbin);      
            else
               SM_tmp=squeeze(SM(i,j,:,:));  
            end

            s1 = s1 + 1;
            h = subplot(varr,varr,s1);
            set(h,'FontSize',13);

            img = squeeze(SM_tmp);
            imagesc([0 30],[-0.2 1],img') % On the frequency axis, 0.5 is equivalent with the Nyquist rate (Fs/2).
            im=gca;
            im.FontSize=14;
%             im.FontWeight='bold';
            colormap jet
            set(h,'YDir','normal')
            caxis([0 1])
            grid on

            if(i==varr && j==ceil(varr/2))
                xlabel('Normalized Frequency','Fontsize',16, 'Fontweight', 'bold')
            end
            if(i==ceil(varr/2) && j==1)
                ylabel('Time (sample)','Fontsize',16, 'Fontweight', 'bold')
            end
            title(['var' num2str(i) ' <--- var' num2str(j)],'Fontsize',14, 'Fontweight', 'bold')
        end
    end
    h = colorbar;
    set(h, 'Position', [.92 .11 .03 .8150],'Fontsize',16, 'Fontweight', 'bold', 'LineWidth', 1)
end
