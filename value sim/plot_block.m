function plot_block(input,tlabel,rule_labels)
% % plot_block %
%PURPOSE:   Plot choice behavior along a block
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   input:  Structure generated by choice_block().
%   tlabel: Text string that will be put on top of the figure
%   rule_labels:    Text strings for the legend in figure

%%
if numel(input)>1  %more than 1 data set, plot mean+-sem
    n=input{1}.n;
    numblockType=input{1}.numblockType;
    blockType=input{1}.blockType;
    for j=1:numel(input)
        temp_probl(:,:,j)=input{j}.probl;
        temp_probr(:,:,j)=input{j}.probr;
        temp_probreward(:,:,j)=input{j}.probreward;
    end
    probl=nanmean(temp_probl,3);
    probr=nanmean(temp_probr,3);
    probreward=nanmean(temp_probreward,3);
    probl_sem=nanstd(temp_probl,[],3)./sqrt(numel(input));
    probr_sem=nanstd(temp_probr,[],3)./sqrt(numel(input));
    probreward_sem=nanstd(temp_probreward,[],3)./sqrt(numel(input));
else        %plot the 1 data set
    n=input.n;
    probl=input.probl;
    probr=input.probr;
    probreward=input.probreward;
    numblockType=input.numblockType;
    blockType=input.blockType;
end

colors = cbrewer('qual', 'Set1', numblockType);

figure;
subplot(1,2,1); hold on;
legstring=[]; 
for j=1:numblockType
    plot(n,probl(:,j),'.-','MarkerSize',30,'Linewidth',3,'Color',colors(j,:));
    legstring{j}=['->' rule_labels{blockType(j)}];
end
if numel(input)>1
    for j=1:numblockType
        for l=1:numel(n)
            plot(n(l)*[1 1],probl(l,j)+probl_sem(l,j)*[-1 1],'-','LineWidth',3,'Color',colors(j,:));
        end
    end
end
legend(legstring);
ylabel('P_{left}');
xlabel('Trials from block switch');
axis([n(1) n(end) 0 1]);
title(tlabel);

subplot(1,2,2); hold on;
for j=1:numblockType
    plot(n,probr(:,j),'.-','MarkerSize',30,'Linewidth',3,'Color',colors(j,:));
end
if numel(input)>1
    for j=1:numblockType
        for l=1:numel(n)
            plot(n(l)*[1 1],probr(l,j)+probr_sem(l,j)*[-1 1],'-','LineWidth',3,'Color',colors(j,:));
        end
    end
end
legend(legstring);
ylabel('P_{right}');
xlabel('Trials from block switch');
axis([n(1) n(end) 0 1]);

print(gcf,'-dpng','blocks');    %png format
saveas(gcf, 'blocks', 'fig');

%% plot prob of reward

figure;

subplot(1,2,1); hold on;
for j=1:numblockType
    plot(n,probreward(:,j),'.-','MarkerSize',30,'Linewidth',3,'Color',colors(j,:));
end
if numel(input)>1
    for j=1:numblockType
        for l=1:numel(n)
            plot(n(l)*[1 1],probreward(l,j)+probreward_sem(l,j)*[-1 1],'-','LineWidth',3,'Color',colors(j,:));
        end
    end
end
legend(legstring);
ylabel('P_{reward}');
xlabel('Trials from block switch');
axis([n(1) n(end) 0 1]);

print(gcf,'-dpng','blocks_reward');    %png format
saveas(gcf, 'blocks_reward', 'fig');

end


    