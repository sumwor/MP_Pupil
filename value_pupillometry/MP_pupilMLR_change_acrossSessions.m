function MP_pupilMLR_change_acrossSessions(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
setup_figprop
set(0, 'DefaultFigureRenderer', 'painters');
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);

% load the first file to get some parameters
all_coeff_future = [];
all_pval_future = [];


all_coeff_iti_1 = [];
all_pval_iti_1 = [];
all_coeff_iti_2 = [];
all_pval_iti_2 = [];
all_coeff_iti_3 = [];
all_pval_iti_3 = [];

for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_change.mat']);
    saveRegName_ITI = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);

    if exist(saveRegName)
        load(saveRegName)
        load(saveRegName_ITI)
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
        
        % load choice and reward regression
%          reg_all.regr_time = reg_cr_change.regr_time;
%         reg_all.numPredictor = reg_cr_change.numPredictor;
%         reg_all.nback = reg_cr_change.nback;
%         reg_all.interaction = reg_cr_change.interaction;
%           all_coeff = cat(3,all_coeff, reg_cr_change.coeff);
%         all_pval = cat(3,all_pval, reg_cr_change.pval);
        
    
    % load the MLR with C(n+1)

        all_coeff_future = cat(3,all_coeff_future, reg_cr_future_change.coeff);
        all_pval_future = cat(3, all_pval_future, reg_cr_future_change.pval);
        
        reg_all.regr_time = reg_cr_future_change.regr_time;
        reg_all.numPredictor = reg_cr_future_change.numPredictor;
        reg_all.nback = reg_cr_future_change.nback;
        reg_all.interaction = reg_cr_future_change.interaction;
     % load the ITI regression (n+1 and n)
        
       
%         all_coeff_iti1 = cat(3,all_coeff_iti1, reg_cr1_change.coeff);
%         all_pval_iti1 = cat(3,all_pval_iti1, reg_cr1_change.pval);
%         
%         all_coeff_iti2 = cat(3,all_coeff_iti2, reg_cr2_change.coeff);
%         all_pval_iti2 = cat(3,all_pval_iti2, reg_cr2_change.pval);
%        
%         all_coeff_iti3 = cat(3,all_coeff_iti3, reg_cr3_change.coeff);
%         all_pval_iti3 = cat(3,all_pval_iti3, reg_cr3_change.pval);
%         
         % load the ITI regression (n-1 and n)
        
        all_coeff_iti_1 = cat(3,all_coeff_iti_1, reg_cr1_change_1.coeff);
        all_pval_iti_1 = cat(3,all_pval_iti_1, reg_cr1_change_1.pval);
        
        all_coeff_iti_2 = cat(3,all_coeff_iti_2, reg_cr2_change_1.coeff);
        all_pval_iti_2 = cat(3,all_pval_iti_2, reg_cr2_change_1.pval);
       
        all_coeff_iti_3 = cat(3,all_coeff_iti_3, reg_cr3_change_1.coeff);
        all_pval_iti_3 = cat(3,all_pval_iti_3, reg_cr3_change_1.pval);

        %%  load control regression for choice and reward regression
        %no control for now
%         
%         if ~exist('all_coeff_ctrl') && exist('reg_cr_change_ctrl')
%             % initialize the control matrix
%             fieldsName = fieldnames(reg_cr_change_ctrl);
%             for tt = 1:length(fieldsName)
%                 all_coeff_ctrl.(fieldsName{tt}) = [];
%                 all_pval_ctrl.(fieldsName{tt}) = [];
%             end
%             
%         end
%         if exist('reg_cr_change_ctrl')
%         % load the regression results
%             fieldsName = fieldnames(reg_cr_change_ctrl);
%             for uu = 1:length(fieldsName)
%                 all_coeff_ctrl.(fieldsName{uu}) = cat(3, all_coeff_ctrl.(fieldsName{uu}), reg_cr_change_ctrl.(fieldsName{uu}).coeff);
%                 all_pval_ctrl.(fieldsName{uu}) = cat(3, all_pval_ctrl.(fieldsName{uu}), reg_cr_change_ctrl.(fieldsName{uu}).pval);
%             end
%             
%         end   
%       
%         % load future control
%         if ~exist('all_coeff_future_ctrl')  && exist('reg_cr_future_change_ctrl')
%             % initialize the control matrix
%             fieldsName_future = fieldnames(reg_cr_future_change_ctrl);
%             for tt = 1:length(fieldsName_future)
%                 all_coeff_future_ctrl.(fieldsName_future{tt}) = [];
%                 all_pval_future_ctrl.(fieldsName_future{tt}) = [];
%             end
%         end
%         % load the regression results
%         
%         if exist('reg_cr_future_change_ctrl')
%             fieldsName_future = fieldnames(reg_cr_future_change_ctrl);
%             for uu = 1:length(fieldsName_future)
%                 all_coeff_future_ctrl.(fieldsName_future{uu}) = cat(3, all_coeff_future_ctrl.(fieldsName_future{uu}), reg_cr_future_change_ctrl.(fieldsName_future{uu}).coeff);
%                 all_pval_future_ctrl.(fieldsName_future{uu}) = cat(3, all_pval_future_ctrl.(fieldsName_future{uu}), reg_cr_future_change_ctrl.(fieldsName_future{uu}).pval);
%             end
%         end
    end
end
        



%% original linear regression
% 1. use bootstrap to get the average and 95% CI for each factor, plot the
% bar plot

%% other things can be done for pupil response:
% correlation: pupil response - latent variable
%                                                     - response time

% reg_all.coeff= all_coeff;
% 
% % use bootstrp to get coefficient
% reg_all = getBootstrp(reg_all, 0, 0.05);
% 

% go to the save path
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

% %figure 1
% figure;
% xAxis = 1:10;
% bar(xAxis(2:10),reg_all.coeff_bootave(2:10),'FaceColor',[0.7,0.7,0.7])                
% 
% hold on
% erneg = reg_all.coeff_bootave(2:10)-reg_all.bootlow(2:10);
% erpos = reg_all.boothigh(2:10)-reg_all.coeff_bootave(2:10);
% er = errorbar(xAxis(2:10),reg_all.coeff_bootave(2:10),erneg,erpos);    
% 
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'xticklabel',{'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'})
% hold off
% pvalThresh = NaN;
% xtickangle(45)
% ylabel('Coefficients (a.u.)');
% title('Coefficient for pupil change - choice and reward');
% pvalThresh = NaN;
% xtitle = 'Time from cue (s)';
% tlabel={'c(n)','c(n-1)','c(n-2)','r(n)','r(n-1)','r(n-2)','c(n)xr(n)','c(n-1)xr(n-1)','c(n-2)xr(n-2)'};
% MP_plot_regrcoef_pupil(reg_all,pvalThresh,tlabel,xtitle);
% 
% print(gcf,'-dpng','MLR-change_choiceoutcome_averageSession');    %png format
% saveas(gcf, 'MLR-change_choiceoutcome_averageSession', 'fig');
% saveas(gcf, 'MLR-change_choiceoutcome_averageSession','svg');
% 
% % plot the figure as number of session that is significant
% reg_sig.coeff = all_coeff;
% reg_sig.pval = all_pval;
% reg_sig.regr_time = reg_all.regr_time;
% reg_sig.numPredictor = reg_all.numPredictor;
% reg_sig.nback = reg_all.nback;
% reg_sig.interaction = reg_all.interaction;
% reg_sig.pvalThresh= 0.01;
% 
% if exist('all_coeff_ctrl')
%     reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01, 0.01);
% 
% % MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
%     MP_plot_regr(reg_sig,reg_pval_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
% else
%     MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
% end
% 
% 
% print(gcf,'-dpng','MLR_change-choiceoutcome_sigSession');    %png format
% saveas(gcf, 'MLR_change-choiceoutcome_sigSession', 'fig');
% saveas(gcf, 'MLR-change_choiceoutcome_sigSession','svg');
% 


%% combine sessions from same animal together
% for ii = 1:length(animalList)
%     newanimalList{ii} = animalList{ii}(1:3);
% end
% newanimalList = unique(newanimalList);
% 
% newsubject_mask = subject_mask;
% for ii=1:length(animalList)
%     if ii == 4
%         newsubject_mask(subject_mask==ii) = 3;
%     elseif ii == 5 || ii == 6 || ii == 7
%         newsubject_mask(subject_mask == ii) = 4;
%     elseif ii == 8
%         newsubject_mask(subject_mask == ii) = 5;
%     end
% end

% %% plot the combined animals
% posSigCn = zeros(size(reg_sig.coeff,1), length(newanimalList));
% negSigCn = zeros(size(reg_sig.coeff,1), length(newanimalList));
% for tt = 1:length(newanimalList)
%     if sum(newsubject_mask == tt) > 0   % if pupil data for certain animal exists
%         reg_sub = reg_sig;
%         reg_sub.coeff = reg_sig.coeff(:,:,newsubject_mask == tt);
%         reg_sub.pval = reg_sig.pval(:,:, newsubject_mask == tt);
%         posSigCn(:,tt) = sum((reg_sub.pval(:,2,:)<reg_sub.pvalThresh&reg_sub.coeff(:,2,:)>0),3)/sum(newsubject_mask == tt);
%         negSigCn(:,tt) = sum((reg_sub.pval(:,2,:)<reg_sub.pvalThresh&reg_sub.coeff(:,2,:)<0),3)/sum(newsubject_mask == tt);
%     end
% end
% 
% figure;
% for uu = 1:length(newanimalList)
%     subplot(5,1,uu)
%     %plot(reg_sig.regr_time, posSig(:,1),'r');
%     H1=area(reg_sig.regr_time,posSigCn(:,uu));
%     set(H1(1),'FaceColor',[1 0.5 0],'EdgeColor',[1,1,0]);
%     hold on;
%     %plot(reg_sig.regr_time, -negSig(:,1));
%     H2=area(reg_sig.regr_time,-negSigCn(:,uu));
%     set(H2(1),'FaceColor',[0 1 0],'EdgeColor',[0,1,0]);
%     hold on;
%     plot([reg_sig.regr_time(1) reg_sig.regr_time(end)], [0 0],'k','LineWidth',1.5);
%     hold on;
%     plot([0 0],[-1 1],'k','LineWidth', 1.5);
%     set(gca,'box','off');
%     set(gca,'XColor','none','YColor','none')
%     pos = get(gca, 'Position');
%     pos(4) = 0.18;
%     set(gca, 'Position', pos)
% end
% print(gcf,'-dpng','MLR_change-posneg5_sigSession');    %png format
% saveas(gcf, 'MLR_change-posneg5_sigSession', 'fig');
% saveas(gcf, 'MLR-change_posneg5_sigSession','svg');
% 
%%  linear regression with C(n+1)
reg_cr_future_change.coeff= all_coeff_future;

% use bootstrp to get coefficient
reg_cr_future_change = getBootstrp(reg_cr_future_change, 0, 0.05);

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);
% 
% figure;
% xAxis = 1:8;
% bar(xAxis(2:8),reg_cr_future.coeff_bootave(2:8),'FaceColor',[0.7,0.7,0.7])                
% 
% hold on
% erneg = reg_cr_future.coeff_bootave(2:8)-reg_cr_future.bootlow(2:8);
% erpos =reg_cr_future.boothigh(2:8)-reg_cr_future.coeff_bootave(2:8);
% er = errorbar(xAxis(2:8),reg_cr_future.coeff_bootave(2:8),erneg,erpos);    
% 
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'xticklabel',{'C(n)','C(n+1)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)'})
% hold off
% 
% xtickangle(45)
% ylabel('Coefficients (a.u.)');
% title('Coefficient for pupil change - choice and reward');
xtitle='Time from cue (s)';
tlabel={'c(n+1)','c(n)','c(n-1)','c(n-2)','r(n+1)','r(n)', 'r(n-1)','r(n-2)',...
        'c(n+1)*r(n+1)','c(n)*r(n)','c(n-1)*r(n-1)','c(n-2)*r(n-2)','Reward Rate','Cumulative Reward'};
pvalThresh=0.01;
MP_plot_regrcoef_pupil(reg_cr_future_change,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-change-choiceoutcome_future_averageSession');    %png format
saveas(gcf, 'MLR-change-choiceoutcome_future_averageSession', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome_future_averageSession','svg');

% plot the figure as number of session that is significant
reg_sig.coeff = all_coeff_future;
reg_sig.pval = all_pval_future;
reg_sig.regr_time = reg_cr_future_change.regr_time;
reg_sig.numPredictor = reg_cr_future_change.numPredictor;
reg_sig.nback = reg_cr_future_change.nback;
reg_sig.interaction = reg_cr_future_change.interaction;
reg_sig.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);

if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR-change-choiceoutcome_future_sigSession');    %png format
saveas(gcf, 'MLR-change-choiceoutcome_future_sigSession', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome_future_sigSession','svg');

%% separate the choice coefficient into pos/neg 

% tlabel={'C(n)-pos','C(n)-neg','C(n)-|pos-neg|'};
% MP_plot_regr_PN(reg_sig,reg_sig.pvalThresh,tlabel,xtitle);
% print(gcf,'-dpng','MLR_change-choiceoutcome_sigSession_abs_pos-neg');    %png format
% saveas(gcf, 'MLR_change-choiceoutcome_sigSession_abs_pos-neg', 'fig');
% saveas(gcf, 'MLR-change_choiceoutcome_sigSession_abs_pos-neg','svg');
% choice effect should be more consistent within subject, not so for motor
% effects?


% % loop through all animals, plot then separatly
% posSigCn = zeros(size(reg_sig.coeff,1), length(animalList));
% negSigCn = zeros(size(reg_sig.coeff,1), length(animalList));
% posSigCn_1 = zeros(size(reg_sig.coeff,1), length(animalList));
% negSigCn_1 = zeros(size(reg_sig.coeff,1), length(animalList));
% for tt = 1:length(animalList)
%     tlabel={'C(n)-pos','C(n)-neg','C(n)-|pos-neg|'};
%     xtitle = ['Time from cue (s)',animalList{tt}] ;
%     if sum(subject_mask == tt) > 0   % if pupil data for certain animal exists
%         reg_sub = reg_sig;
%         reg_sub.coeff = reg_sig.coeff(:,:,subject_mask == tt);
%         reg_sub.pval = reg_sig.pval(:,:, subject_mask == tt);
%         posSigCn(:,tt) = sum((reg_sub.pval(:,3,:)<reg_sub.pvalThresh&reg_sub.coeff(:,3,:)>0),3)/sum(subject_mask == tt);
%         negSigCn(:,tt) = sum((reg_sub.pval(:,3,:)<reg_sub.pvalThresh&reg_sub.coeff(:,3,:)<0),3)/sum(subject_mask == tt);
%         posSigCn_1(:,tt) = sum((reg_sub.pval(:,4,:)<reg_sub.pvalThresh&reg_sub.coeff(:,4,:)>0),3)/sum(subject_mask == tt);
%         negSigCn_1(:,tt) = sum((reg_sub.pval(:,4,:)<reg_sub.pvalThresh&reg_sub.coeff(:,4,:)<0),3)/sum(subject_mask == tt);
%         posSigRn(:,tt) = sum((reg_sub.pval(:,7,:)<reg_sub.pvalThresh&reg_sub.coeff(:,7,:)>0),3)/sum(subject_mask == tt);
%         negSigRn(:,tt) = sum((reg_sub.pval(:,7,:)<reg_sub.pvalThresh&reg_sub.coeff(:,7,:)<0),3)/sum(subject_mask == tt);
%    
%         %MP_plot_regr_PN(reg_sub,reg_sub.pvalThresh,tlabel,xtitle);
%         %print(gcf,'-dpng',['MLR-change-choiceoutcome_sigSession_abs_pos-neg',animalList{tt}]);    %png format
%         %saveas(gcf, ['MLR-change-choiceoutcome_sigSession_abs_pos-neg',animalList{tt}], 'fig');
% 
%     end
% end
% 
% % Cn
% figure;
% for uu = 1:length(animalList)
%     subplot(length(animalList),1,uu)
%     %plot(reg_sig.regr_time, posSig(:,1),'r');
%     H1=area(reg_sig.regr_time,posSigCn(:,uu));
%     set(H1(1),'FaceColor',[1 0.5 0],'EdgeColor',[1,1,0]);
%     hold on;
%     %plot(reg_sig.regr_time, -negSig(:,1));
%     H2=area(reg_sig.regr_time,-negSigCn(:,uu));
%     set(H2(1),'FaceColor',[0 1 0],'EdgeColor',[0,1,0]);
%     hold on;
%     plot([reg_sig.regr_time(1) reg_sig.regr_time(end)], [0 0],'k','LineWidth',1.5);
%     hold on;
%     plot([0 0],[-1 1],'k','LineWidth', 1.5);
%     set(gca,'box','off');
%     set(gca,'XColor','none','YColor','none')
%     pos = get(gca, 'Position');
%     pos(4) = 0.14;
%     set(gca, 'Position', pos)
% end
% sgtitle('Choice coefficient in subjects');
% print(gcf,'-dpng','MLR_change-posnegCn_sigSession');    %png format
% saveas(gcf, 'MLR_change-posnegCn_sigSession', 'fig');
% saveas(gcf, 'MLR-change_posnegCn_sigSession','svg');
% 
% load the data into one vector
%[x,y] = size(posSigCn);
%CnPosRe = reshape(posSigCn,x*y,1); 
%CnNegRe = reshape(negSigCn,x*y,1); 


%pPosCn = vartestn(CnPosRe,AnimalMask,'TestType','LeveneAbsolute');
%pNegCn = vartestn(CnNegRe,AnimalMask,'TestType','LeveneAbsolute');
% Cn-1
% figure;
% for uu = 1:length(animalList)
%     subplot(length(animalList),1,uu)
%     %plot(reg_sig.regr_time, posSig(:,1),'r');
%     H1=area(reg_sig.regr_time,posSigCn_1(:,uu));
%     set(H1(1),'FaceColor',[1 0.5 0],'EdgeColor',[1,1,0]);
%     hold on;
%     %plot(reg_sig.regr_time, -negSig(:,1));
%     H2=area(reg_sig.regr_time,-negSigCn_1(:,uu));
%     set(H2(1),'FaceColor',[0 1 0],'EdgeColor',[0,1,0]);
%     hold on;
%     plot([reg_sig.regr_time(1) reg_sig.regr_time(end)], [0 0],'k','LineWidth',1.5);
%     hold on;
%     plot([0 0],[-1 1],'k','LineWidth', 1.5);
%     set(gca,'box','off');
%     set(gca,'XColor','none','YColor','none')
%     pos = get(gca, 'Position');
%     pos(4) = 0.14;
%     set(gca, 'Position', pos)
% end
% sgtitle('Past choice coefficient in subjects');
% print(gcf,'-dpng','MLR_change-posnegCn_1_sigSession');    %png format
% saveas(gcf, 'MLR_change-posnegCn_1_sigSession', 'fig');
% saveas(gcf, 'MLR-change_posnegCn_1_sigSession','svg');
% 
% %Rn
% figure;
% 
% for uu = 1:length(animalList)
%     subplot(length(animalList),1,uu)
%     %plot(reg_sig.regr_time, posSig(:,1),'r');
%     H1=area(reg_sig.regr_time,posSigRn(:,uu));
%     set(H1(1),'FaceColor',[1 0.5 0],'EdgeColor',[1,1,0]);
%     hold on;
%     %plot(reg_sig.regr_time, -negSig(:,1));
%     H2=area(reg_sig.regr_time,-negSigRn(:,uu));
%     set(H2(1),'FaceColor',[0 1 0],'EdgeColor',[0,1,0]);
%     hold on;
%     plot([reg_sig.regr_time(1) reg_sig.regr_time(end)], [0 0],'k','LineWidth',1.5);
%     hold on;
%     plot([0 0],[-1 1],'k','LineWidth', 1.5);
%     set(gca,'box','off');
%     set(gca,'XColor','none','YColor','none')
%     pos = get(gca, 'Position');
%     pos(4) = 0.14;
%     set(gca, 'Position', pos)
% end
% sgtitle('Outcome coefficient in subjects');
% print(gcf,'-dpng','MLR_change-posnegRn_sigSession');    %png format
% saveas(gcf, 'MLR_change-posnegRn_sigSession', 'fig');
% saveas(gcf, 'MLR-change_posnegRn_sigSession','svg');
% 
% % use m box test
% % group1: choice; group2: Rn (positive)
% TestMatCnPos = [ones(1,length(animalList));posSigCn];
% TestMatRnPos = [ones(1,length(animalList))*2;posSigRn];
% TestMatPos = [TestMatCnPos,TestMatRnPos]';
% MBoxtest(TestMatPos',0.05)
% 
% %% test for same variance for each time point
% posP = zeros(1,size(posSigRn,1));
% group = [ones(1,size(posSigRn,2)),ones(1,size(posSigRn,2))*2];
% for tt = 1:size(posSigRn,1)
%     posP(tt) = vartestn([posSigCn(tt,:), posSigRn(tt,:)]',group','TestType','LeveneAbsolute','Display','off');
% end
% figure;plot(reg_all.regr_time,posP)
% hold on; plot([reg_all.regr_time(1) reg_all.regr_time(end)],[0.05 0.05])
% title('Positive coefficients')
% ylabel('p value');
% xlabel('Time from cue(s)');
% print(gcf,'-dpng','MLR_change-pos_vartest');    %png format
% saveas(gcf, 'MLR_change-pos_vartest', 'fig');
% saveas(gcf, 'MLR-change_pos_vartest','svg');
% 
% 
% negP = zeros(1,size(posSigRn,1));
% group = [ones(1,size(posSigRn,2)),ones(1,size(posSigRn,2))*2];
% for tt = 1:size(posSigRn,1)
%     negP(tt) = vartestn([negSigCn(tt,:), negSigRn(tt,:)]',group','TestType','LeveneAbsolute','Display','off');
% end
% figure;plot(reg_all.regr_time,negP)
% hold on; plot([reg_all.regr_time(1) reg_all.regr_time(end)],[0.05 0.05])
% title('Negative coefficients')
% ylabel('p value');
% xlabel('Time from cue(s)');
% print(gcf,'-dpng','MLR_change-neg_vartest');    %png format
% saveas(gcf, 'MLR_change-neg_vartest', 'fig');
% saveas(gcf, 'MLR-change_neg_vartest','svg');
% 
% %% get average coefficient, test for equal variance
% CnCoeff = zeros(size(reg_sig.coeff,1), length(animalList));
% RnCoeff = zeros(size(reg_sig.coeff,1), length(animalList));
% 
% for tt = 1:length(animalList)
%     tlabel={'C(n)-pos','C(n)-neg','C(n)-|pos-neg|'};
%     xtitle = ['Time from cue (s)',animalList{tt}] ;
%     if sum(subject_mask == tt) > 0   % if pupil data for certain animal exists
%         reg_sub = reg_sig;
%         reg_sub.coeff = reg_sig.coeff(:,:,subject_mask == tt);
%         reg_sub.pval = reg_sig.pval(:,:, subject_mask == tt);
%         CnCoeff(:,tt) = sum(reg_sub.coeff(:,3,:),3)/sum(subject_mask == tt);
%         RnCoeff(:,tt) = sum(reg_sub.coeff(:,7,:),3)/sum(subject_mask == tt);
%        
%         %MP_plot_regr_PN(reg_sub,reg_sub.pvalThresh,tlabel,xtitle);
%         %print(gcf,'-dpng',['MLR-change-choiceoutcome_sigSession_abs_pos-neg',animalList{tt}]);    %png format
%         %saveas(gcf, ['MLR-change-choiceoutcome_sigSession_abs_pos-neg',animalList{tt}], 'fig');
% 
%     end
% end
% 
% 
% PCheck = zeros(1,size(CnCoeff,1));
% group = [ones(1,size(CnCoeff,2)),ones(1,size(CnCoeff,2))*2];
% for tt = 1:size(CnCoeff,1)
%     PCheck(tt) = vartestn([CnCoeff(tt,:), RnCoeff(tt,:)]',group','TestType','LeveneAbsolute','Display','off');
% end
% figure;plot(reg_all.regr_time,PCheck)
% hold on; plot([reg_all.regr_time(1) reg_all.regr_time(end)],[0.05 0.05])
% title('Equal variance test for mean coefficient')
% ylabel('p value');
% xlabel('Time from cue(s)');
% print(gcf,'-dpng','MLR_change-mean_vartest');    %png format
% saveas(gcf, 'MLR_change-mean_vartest', 'fig');
% saveas(gcf, 'MLR-change_mean_vartest','svg');
% 
% CnCoeffMedian = zeros(size(reg_sig.coeff,1), length(animalList));
% RnCoeffMedian = zeros(size(reg_sig.coeff,1), length(animalList));
% 
% for tt = 1:length(animalList)
%     tlabel={'C(n)-pos','C(n)-neg','C(n)-|pos-neg|'};
%     xtitle = ['Time from cue (s)',animalList{tt}] ;
%     if sum(subject_mask == tt) > 0   % if pupil data for certain animal exists
%         reg_sub = reg_sig;
%         reg_sub.coeff = reg_sig.coeff(:,:,subject_mask == tt);
%         reg_sub.pval = reg_sig.pval(:,:, subject_mask == tt);
%         CnCoeffMedian(:,tt) = median(reg_sub.coeff(:,3,:),3);
%         RnCoeffMedian(:,tt) = median(reg_sub.coeff(:,7,:),3);
%        
%         %MP_plot_regr_PN(reg_sub,reg_sub.pvalThresh,tlabel,xtitle);
%         %print(gcf,'-dpng',['MLR-change-choiceoutcome_sigSession_abs_pos-neg',animalList{tt}]);    %png format
%         %saveas(gcf, ['MLR-change-choiceoutcome_sigSession_abs_pos-neg',animalList{tt}], 'fig');
% 
%     end
% end
% 
% 
% PCheckMedian = zeros(1,size(CnCoeffMedian,1));
% group = [ones(1,size(CnCoeffMedian,2)),ones(1,size(CnCoeffMedian,2))*2];
% for tt = 1:size(CnCoeff,1)
%     PCheckMedian(tt) = vartestn([CnCoeffMedian(tt,:), RnCoeffMedian(tt,:)]',group','TestType','LeveneAbsolute','Display','off');
% end
% figure;plot(reg_all.regr_time,PCheckMedian, 'k')
% hold on; plot([reg_all.regr_time(1) reg_all.regr_time(end)],[0.05 0.05])
%  set(gca,'box','off');
% title('Equal variance test for median coefficients')
% ylabel('p value');
% xlabel('Time from cue(s)');
% print(gcf,'-dpng','MLR_change-median_vartest');    %png format
% saveas(gcf, 'MLR_change-median_vartest', 'fig');
% saveas(gcf, 'MLR-change_median_vartest','svg');
% 
% % plot the median coefficient for every animal
% figure;
% 
% for uu = 1:length(animalList)
%     subplot(length(animalList),1,uu)
%     plot([5 5], [0 1], 'k', 'LineWidth',2);
%     hold on;
%     %plot(reg_sig.regr_time, posSig(:,1),'r');
%     H1=area(reg_sig.regr_time,CnCoeffMedian(:,uu));
%     set(H1(1),'FaceColor',[1 0.5 0],'EdgeColor',[1,1,0]);
%     hold on;
%     plot([reg_sig.regr_time(1) reg_sig.regr_time(end)], [0 0],'k','LineWidth',1.5);
%     hold on;
%     plot([0 0],[-1 1],'k','LineWidth', 0.5);
%     set(gca,'box','off');
%     if uu<length(animalList) 
%         set(gca,'XColor','none','YColor','none');
%     else
%         set(gca,'YColor','none');
%     end
%    
%     pos = get(gca, 'Position');
%     pos(4) = 0.14;
%     set(gca, 'Position', pos)
%     
%       
% end
% sgtitle('Median coefficient of current choice');
% print(gcf,'-dpng','MLR_change-medianCn');    %png format
% saveas(gcf, 'MLR_change-medianCn', 'fig');
% saveas(gcf, 'MLR_change-medianCn','svg');
% 
% figure;
% 
% for uu = 1:length(animalList)
%     subplot(length(animalList),1,uu)
%     plot([5 5], [0 1], 'k', 'LineWidth',2);
%     hold on;
%     %plot(reg_sig.regr_time, posSig(:,1),'r');
%     H1=area(reg_sig.regr_time,RnCoeffMedian(:,uu));
%     set(H1(1),'FaceColor',[0 1 0],'EdgeColor',[1,1,0]);
%     hold on;
%     plot([reg_sig.regr_time(1) reg_sig.regr_time(end)], [0 0],'k','LineWidth',1.5);
%     hold on;
%     plot([0 0],[-1 1],'k','LineWidth', 0.5);
%     set(gca,'box','off');
%      if uu<length(animalList) 
%         set(gca,'XColor','none','YColor','none');
%     else
%         set(gca,'YColor','none');
%     end
%     pos = get(gca, 'Position');
%     pos(4) = 0.14;
%     set(gca, 'Position', pos)
% end
% sgtitle('Median coefficient of current reward');
% print(gcf,'-dpng','MLR_change-medianRn');    %png format
% saveas(gcf, 'MLR_change-medianRn', 'fig');
% saveas(gcf, 'MLR_change-medianRn','svg');
%pPosRn = vartestn(RnPosRe,AnimalMask,'TestType','LeveneAbsolute');
%pNegRn = vartestn(RnNegRe,AnimalMask,'TestType','LeveneAbsolute');


%% plot the coefficient in different ITIs: aligned to the ITI between n+1 and n
% % plot the figure as number of session that is significant
% reg_sig_iti1.coeff = all_coeff_iti1;
% reg_sig_iti1.pval = all_pval_iti1;
% reg_sig_iti1.regr_time = reg_cr1_change.regr_time;
% reg_sig_iti1.numPredictor = reg_cr1_change.numPredictor;
% reg_sig_iti1.nback = reg_cr1_change.nback;
% reg_sig_iti1.interaction = reg_cr1_change.interaction;
% reg_sig_iti1.pvalThresh= 0.01;
% 
% 
% reg_sig_iti2.coeff = all_coeff_iti2;
% reg_sig_iti2.pval = all_pval_iti2;
% reg_sig_iti2.regr_time = reg_cr2_change.regr_time;
% reg_sig_iti2.numPredictor = reg_cr2_change.numPredictor;
% reg_sig_iti2.nback = reg_cr2_change.nback;
% reg_sig_iti2.interaction = reg_cr2_change.interaction;
% reg_sig_iti2.pvalThresh= 0.01;
% 
% reg_sig_iti3.coeff = all_coeff_iti3;
% reg_sig_iti3.pval = all_pval_iti3;
% reg_sig_iti3.regr_time = reg_cr3_change.regr_time;
% reg_sig_iti3.numPredictor = reg_cr3_change.numPredictor;
% reg_sig_iti3.nback = reg_cr3_change.nback;
% reg_sig_iti3.interaction = reg_cr3_change.interaction;
% reg_sig_iti3.pvalThresh= 0.01;
% 
% tlabel={'C(n+1)'};
% 
% xtitle='Time from cue (s)';
% MP_plot_regr_iti(reg_sig_iti1, reg_sig_iti2, reg_sig_iti3, reg_sig.pvalThresh,tlabel,xtitle);
% print(gcf,'-dpng','MLR-change-choiceoutcome_ITI');    %png format
% saveas(gcf, 'MLR-change-choiceoutcome_ITI', 'fig');
% saveas(gcf, 'MLR-change-choiceoutcome_ITI','svg');
% 
% tlabel={'C(n)'};
% 
% MP_plot_regr_iti(reg_sig_iti1, reg_sig_iti2, reg_sig_iti3, reg_sig.pvalThresh,tlabel,xtitle);
% 
% print(gcf,'-dpng','MLR-change_ITI_Cn');    %png format
% saveas(gcf, 'MLR-change_ITI_Cn', 'fig');
% saveas(gcf, 'MLR-change_ITI_Cn','svg');
% 
% tlabel={'R(n-1)'};
% MP_plot_regr_iti(reg_sig_iti1, reg_sig_iti2, reg_sig_iti3, reg_sig.pvalThresh,tlabel,xtitle);
% print(gcf,'-dpng','MLR-change_ITI-R-1');    %png format
% saveas(gcf, 'MLR-change_ITI-R-1', 'fig');
% saveas(gcf, 'MLR-change_ITI-R-1','svg');

reg_sig_iti_1.coeff = all_coeff_iti_1;
reg_sig_iti_1.pval = all_pval_iti_1;
reg_sig_iti_1.regr_time = reg_cr1_change_1.regr_time;
reg_sig_iti_1.numPredictor = reg_cr1_change_1.numPredictor;
reg_sig_iti_1.nback = reg_cr1_change_1.nback;
reg_sig_iti_1.interaction = reg_cr1_change_1.interaction;
reg_sig_iti_1.pvalThresh= 0.01;


reg_sig_iti_2.coeff = all_coeff_iti_2;
reg_sig_iti_2.pval = all_pval_iti_2;
reg_sig_iti_2.regr_time = reg_cr2_change_1.regr_time;
reg_sig_iti_2.numPredictor = reg_cr2_change_1.numPredictor;
reg_sig_iti_2.nback = reg_cr2_change_1.nback;
reg_sig_iti_2.interaction = reg_cr2_change_1.interaction;
reg_sig_iti_2.pvalThresh= 0.01;

reg_sig_iti_3.coeff = all_coeff_iti_3;
reg_sig_iti_3.pval = all_pval_iti_3;
reg_sig_iti_3.regr_time = reg_cr3_change_1.regr_time;
reg_sig_iti_3.numPredictor = reg_cr3_change_1.numPredictor;
reg_sig_iti_3.nback = reg_cr3_change_1.nback;
reg_sig_iti_3.interaction = reg_cr3_change_1.interaction;
reg_sig_iti_3.pvalThresh= 0.01;

tlabel={'C(n-1)'};
MP_plot_regr_iti(reg_sig_iti_1, reg_sig_iti_2, reg_sig_iti_3, reg_sig.pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-change-choiceoutcome_ITI_1');    %png format
saveas(gcf, 'MLR-change-choiceoutcome_ITI_1', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome_ITI_1','svg');
%%
close all

end