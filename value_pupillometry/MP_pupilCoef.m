function MP_pupilCoef(dataIndex, savefigpath)

% plot the coefficient of chosen q and R to determine whether the session
% encodes RPE or updated chosen Q
nFiles = size(dataIndex,1);


% try time = 2.45 second first
coeff_chosenQ = [];
coeff_R = [];
coeff_chosenQP = [];
coeff_RP = [];
subject_mask = [];
animalList = unique(dataIndex.Animal);
t_coeff = 2.45;
for ii = 1:nFiles
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    saveRegName_cue = fullfile(fn_beh.folder,'analysis-pupil',[fn_beh.name(1:end-7),'regRL.mat']);
    
    animal = fn_beh.name(1:3);
    
   
    if exist(saveRegName_cue)
        subject_mask = [subject_mask, find(contains(animalList, animal))];
        load(saveRegName_cue)
    % the R is the second factor, chosen q is the 5th factor
        coeff_chosenQ = [coeff_chosenQ;reg_cr.coeff(:,5)'];
        coeff_R =[coeff_R; reg_cr.coeff(:,2)'];
        coeff_chosenQP =[coeff_chosenQP; reg_cr.pval(:,5)'];
        coeff_RP = [coeff_RP,; reg_cr.pval(:,2)'];
    end
        
end


plot1 = coeff_chosenQ(:,(find((reg_cr.regr_time<t_coeff),1,'last')));
pval1 = coeff_chosenQP(:,(find((reg_cr.regr_time<t_coeff),1,'last')));
plot2 = coeff_R(:,(find((reg_cr.regr_time<t_coeff),1,'last')));
pval2 = coeff_RP(:,(find((reg_cr.regr_time<t_coeff),1,'last')));
figure; 
cmap = colormap('prism');
for ii = 1:length(animalList)  % plot different animals in different color
    index_sig = (subject_mask==ii)' & pval1<0.05 & pval2<0.05;
    index_unsig = (subject_mask==ii)' & (pval1>0.05 | pval2>0.05);
    scatter(plot1(index_sig),plot2(index_sig),80,cmap(ii*2,:),'filled');
    hold on;
    scatter(plot1(index_unsig),plot2(index_unsig),80,cmap(ii*2,:));
    hold on;
end
ax = gca;
xlabel('Coefficient of chosen value');
ylabel('Coefficient of R(n)');
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%legend(animalList)
title(['Coefficient of R(n) and chosen Q at t = ',num2str(t_coeff)]);
print(gcf,'-dpng',['coeff_r_chosenq_t=',num2str(t_coeff*100)]);    %png format
saveas(gcf, ['coeff_r_chosenq_t=',num2str(t_coeff*100)], 'fig');
