function MP_plot_regrcoef_RLpupil(input,pvalThresh,tlabel,xtitle)
% % plot_regr %
%PURPOSE:   Plot results from multiple linear regression
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   input:        Structure generated by linear_regr().
%   pvalThresh:   Threshold value to deem whether p-values are significant
%   tlabel:       Text to put as title of the plot.
%   xtitle:       Text to put as the label for x-axis.

%% setup
t=input.regr_time;
dt=nanmean(diff(t));

nCells=numel(input);
for j=1:nCells
    if isfield(input, 'pval')
        pval(:,:,j)=input.pval;
    end
    coeff(:,:,j) = nanmean(input.coeff,3);
end

nPredictor=input.numPredictor;


if (input.interaction == true)
    if nPredictor == 2
        panelv = nPredictor + 1;    %plot extra row for the interaction terms
        nInteraction = 1;
    else
        panelv = nPredictor + 2;    %plot extra row for the interaction terms
        nInteraction = 2;
    end
else
    panelv = nPredictor;
    nInteraction = 0;
end

%% plot results

figure;

for l=1:nPredictor
        currPredictor=1+l; %first +1 because first term is bias
        
        subplot(4,3,currPredictor-1); hold on;
        % patch([t(1) t(end) t(end) t(1)],[0 0 100*pvalThresh 100*pvalThresh],[0.7 0.7 0.7],'EdgeColor','none');
        plot(t,coeff(:,currPredictor,:),'k.-','MarkerSize',30);
        
        xlim([t(1) t(end)]);
        title(tlabel{currPredictor-1});
        xlabel(xtitle);
        ylabel('Coefficients');
        yl = ylim;
        % get significant point
        if isfield(input, 'pval')  % if there is a pvalue
            for ii = 1:length(pval(:,currPredictor,:))
                if pval(ii,currPredictor,:) < pvalThresh
                    hold on;
                    scatter(t(ii),ceil(yl(2)/2),25,'filled','black');
                end
            end
           
        else
            % this is average across session, obtained a CI with bootstrap
            % plot the shaded errors
            hold on;
            gray=[0.7 0.7 0.7];
            errorshade(t,input.bootlow(:,currPredictor),input.boothigh(:,currPredictor),gray);
            hold on;
            plot(t,coeff(:,currPredictor,:),'k.-','MarkerSize',30);
            
        end
        % plot the vertical line aligned to the cue
        yl = ylim;
        plot([0 0],[-0.5 1],'k','LineWidth',1);
   
end

if nInteraction > 0
    for l = 1:nInteraction
        for k=1:1+nback
            
            currPredictor=1+nPredictor*(1+nback)+(l-1)*(1+nback)+k;
            
            subplot(panelv,1+nback,currPredictor-1); hold on;
            
            plot(t,coeff(:,currPredictor,:),'k.-','MarkerSize',30);
            
            xlim([t(1) t(end)]);
            title(tlabel{currPredictor-1});
            xlabel(xtitle);
            ylabel('Coefficients');
            yl = ylim;
            % get significant point
            if isfield(input, 'pval')
                for ii = 1:length(pval(:,currPredictor,:))
                    if pval(ii,currPredictor,:) < pvalThresh
                        hold on;
                        scatter(t(ii),ceil(yl(2)),25,'filled','black');
                    end
                end
              
            else
                hold on;
                errorshade(t,input.bootlow(:,currPredictor),input.boothigh(:,currPredictor),gray);
                hold on;
                plot(t,coeff(:,currPredictor,:),'k.-','MarkerSize',30);
             
            end
             % plot the vertical line aligned to the cue
        	yl = ylim;
            plot([0 0],[-0.5 0.5],'k','LineWidth',1);
        end
    end
end

end
