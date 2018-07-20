function plottingTF(outputname)
% Plotting PDC-graphs
% Before using it, you should use 'GeneralLinearKalmanFilteringUni.m' to
%   calculate all necessary parameters for AutoRegressive model
% After that you should use 'PDC_calculating.m' to alculate all necessary
%   parameters for PDC

% GeneralLinearKalmanFilteringUni
% PDC_calculating

% clear all;
clearvars -except outputname; 
close all
% Write what data you want to send in PDC { A - Abstractness_TimeSeries
%                                           CmA - CommonAbstractness_TimeSeries
%                                           CmC - CommonConcreteness_TimeSeries
%                                           C - Concreteness_TimeSeries         }

load(outputname)

fh=figure('Name', 'PDC', 'visible', 'off', 'Position', get(0, 'Screensize'));  % 
gh=figure('Name', 'DTF', 'visible', 'off', 'Position', get(0, 'Screensize'));    % 
axis([0 FreqBrainMax 0 1])
% x = linspace(0,FreqBrainMax,size(f,2));
for i=1:m
    for j=1:m
        for timeV=1:size(PDC,1)
            for freqV=1:size(PDC,2)
                PDCplt(timeV, freqV) = PDC{timeV,freqV}(i,j);
                DTFplt(timeV, freqV) = DTF{timeV,freqV}(i,j);
            end
        end
        
%         fh = findobj( 'Type', 'Figure', 'Name', 'PDC' );
        set(0,'CurrentFigure', fh);
        subplot( m,m,((i-1)*m)+j);
        contourf(PDCplt);
        ax = gca;
        ax.XTick = [linspace(0,size(PDCplt,2),6)];
        set(ax,'XTickLabel',linspace(0,FreqBrainMax,6)) 
%         set(ax,'XTickLabel',[]) 
%         xlabel('Frequency, Hz');
        ax.YTick = [linspace(0,size(PDCplt,1),18)];
        set(ax,'YTickLabel',round(linspace(-0.2,1.5,18),1)); 
%         set(ax,'YTickLabel',[]); 
%         ylabel('Time, s');
        
%         gh = findobj( 'Type', 'Figure', 'Name', 'DTF' );
        set(0,'CurrentFigure', gh);
        subplot( m,m,((i-1)*m)+j);
        contourf(DTFplt);
        ax = gca;
        ax.XTick = [linspace(0,size(PDCplt,2),6)];
        set(ax,'XTickLabel',linspace(0,FreqBrainMax,6)) 
%         set(ax,'XTickLabel',[])        
%         xlabel('Frequency, Hz');
        ax.YTick = [linspace(0,size(PDCplt,1),18)];
        set(ax,'YTickLabel',round(linspace(-0.2,1.5,18),1)); 
%         set(ax,'YTickLabel',[]);%         ylabel('Time, s');
        
%             ylabel(['\mid\pi_{',int2str(j),'\rightarrow',num2str(i),'}\mid^2']);
%             title(['Channel y_{',int2str(j),'}\rightarrowy_{',num2str(i),'}']);
        
    end
end
idcs   = strfind(outputname,'\');
newdir = outputname(1:idcs(end-1)-1);
name1=outputname(idcs(end-1)+1:idcs(end)-1);
name2=outputname(idcs(end)+1:end);

figurenamePDC=[newdir, '\PDC', name1,name2,'.png'];
F=getframe(fh);
imwrite(F.cdata, figurenamePDC, 'png')
% set(0,'CurrentFigure', fh);
% print(figurenamePDC,'-dpng')
% colorbar
figurenameDTF=[newdir, '\DTF', name1,name2,'.png'];
G=getframe(gh);
imwrite(G.cdata, figurenameDTF, 'png')
% set(0,'CurrentFigure', gh);
% print(figurenameDTF,'-dpng')
% colorbar

% end
