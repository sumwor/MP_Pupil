% Master file for analyzing pupillometry data for matchingpennies task

%
% To run this code:
% 1) Add all the subfolders to the Path in MATLAB
% 2) Change the variable 'root_path' below
%
% List of unresolved issues:
% 1) Why do the response times have two peaks?


clearvars;
close all;
setup_figprop;

%root_path = '/Volumes/haha/MatchingPennies/pupilData';
root_path = 'E:\labcode\MP_Pupil\pupilData';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Preprocessing: parse each log file and tabulate the data set

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('-----------------------------------------------------------');

% Look for data files and create a database index
logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex = MP_pupil_createBehMatFiles(dataIndex);

% Add information about lesion
% Future - Add information about pupil, imaging, etc.
dataIndex = addIndexPupil(dataIndex);

dataIndex = sortdataIndex(dataIndex);

% Determine if each session fulfill performance criteria (Matching Pennies)
MP_determineBehCriteria(dataIndex);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - single session

% exampleLogName = '19107_phase3_R71NoCue_1904051629.log';
% idx = find(ismember(dataIndex.LogFileName,exampleLogName)==1);

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    savematpath = [fullfile(dataIndex.BehPath{ii},dataIndex.LogFileName{ii}(1:end-4))];
    MP_session(dataIndex.BehPath{ii},dataIndex.LogFileName{ii},savematpath);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - summary
% save_path = fullfile(root_path,'figs_summary');
% bandit_bayesian(dataIndex, save_path)


save_path = fullfile(root_path,'summary','figs_summary');


MP_behaviorPerAnimal(dataIndex,save_path);

MP_behaviorAll(dataIndex(dataIndex.pupil==1,:), save_path);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - Model fitting

model_path = fullfile(root_path,'mat_models');


MP_fittingPerAnimal(dataIndex,model_path);

% save the predicted latent variable into each behavior mat files
MP_saveLatent(dataIndex, model_path);

%% fit for changing parameters
MP_fittingDrift(dataIndex,model_path);
MP_plotDrift(dataIndex,model_path);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - Model comparison

model_path = fullfile(root_path,'mat_models');
MP_compareModels(model_path);   % model comparison

%MP_confusionMat(model_path);
MP_choicekernel_switch(dataIndex(dataIndex.pupil==1,:),model_path)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  pupil analysis

% get the raw csv data, preprocess it

createPupilFiles(dataIndex(dataIndex.pupil==1,:));

%% save path
model_path = fullfile(root_path,'mat_models');
save_path_pupil = fullfile(root_path,'summary','figs_summary_pupil');

%% pupil simple plots
MP_pupilSimpleplots(dataIndex(dataIndex.pupil==1,:));

%% tonic activity, no deterministic results
% MP_tonic(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

% use latent variables to predict 

%% calculate pupil change/pupil response
MP_pupilChange(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

%% fit the softmax to pupil response
MP_pupilResp_latent(dataIndex(dataIndex.pupil==1,:),save_path_pupil,model_path);

%% Linear regression with choice and reward
% running regression (choice and reward) on individual sessions
MP_pupilMLR(dataIndex(dataIndex.pupil==1,:));

% a regression with C(n+1), included in the MP_pupilMLR
%MP_pupilMLR_future(dataIndex(dataIndex.pupil==1,:));
% checking the effect of long term reward
% MP_pupilReward(dataIndex(dataIndex.pupil==1,:))
% MP_pupilReward_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

% regression with pupil position? (saccade)
%take individual regressions plot them together
% MP_pupilMLR_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

MP_pupilMLR_change_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
% multilinear regression for animals
% try average across sessions

%% linear regression with latent variables
% FQ_RPE algorithm5
%MP_pupilRL_MLR(dataIndex(dataIndex.pupil==1,:));
% MP_pupilRL_MLR_all(dataIndex,save_path_pupil);
%MP_pupilRL_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

%MP_pupilRL_change_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

% FQ_RPE_CK algorithm
MP_pupilRL_MLR_withCK(dataIndex(dataIndex.pupil==1,:));
MP_pupilRL_change_acrossSessions_withCK(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
%MP_pupilRL_acrossSessions_withCK(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

% choice vs value
MP_pupilRL_delta_MLR_withCK(dataIndex(dataIndex.pupil==1,:));
MP_pupilRL_acrossSessions_delta(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
%% linear regression to get the RPE or Q
% MP_pupilRL_RPE(dataIndex(dataIndex.pupil==1,:));
% MP_pupilRLRPE_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
% MP_pupilRLRPE_change_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
% 
% scatter plot of the coefficients of chosen Q and R
% MP_pupilCoef(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
% MP_pupilCoef_change(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

% choice autocorrelation
MP_pupilRL_RPE_CK(dataIndex(dataIndex.pupil==1,:));
% MP_pupilRLRPE_acrossSessions_CK(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
MP_pupilRLRPE_change_acrossSessions_CK(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

MP_pupilRL_CKE_CK(dataIndex(dataIndex.pupil==1,:));
MP_pupilRLCKE_change_acrossSessions_CK(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

MP_pupilRL_RPE_rn_CK(dataIndex(dataIndex.pupil==1,:));
MP_pupilRLRPE_rn_change_acrossSessions_CK(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

%% behavior relate to high/low pupil diameter
MP_arousal_beh(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
%% run linear regression with trials separated into high/low pupil diameter
MP_pupilMLR_arousal(dataIndex(dataIndex.pupil==1,:));
MP_pupilMLR_arousal_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

MP_pupil_AS_arousal(dataIndex(dataIndex.pupil==1,:));
MP_pupil_AS_arousal_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

MP_pupil_RPE_arousal(dataIndex(dataIndex.pupil==1,:));
% MP_pupilRLRPE_acrossSessions_CK(dataIndex(dataIndex.pupil==1,:), save_path_pupil);
MP_pupil_RPE_arousal_acrossSessions(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Simulation - using algorithm and actual choices and outcomes from session

saveSimPath = fullfile(root_path,'summary','sim');
if ~exist(saveSimPath)
    mkdir(saveSimPath)
end

exampleLogName = '19107_phase3_R71NoCue_1904051629.log';
exampleAnimal = '19107';
model_path = fullfile(root_path,'mat_models/');
model_type = 'FQ_RPE';

savebehfigpath = bandit_predictSession(dataIndex(normalSubset & reversalSubset,:),exampleLogName,exampleAnimal,model_path,model_type);

print(gcf,'-dpng',fullfile(savebehfigpath,['session_' model_type]));
saveas(gcf, fullfile(savebehfigpath,['session_' model_type]), 'fig');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Simulation - using algorithm and task design

n_stim=10000;    % number of trials to simulate


% using median value from fit with animal data
model_path = fullfile(root_path,'mat_models');
AllModelfiles = dir(fullfile(model_path,'*_model.mat'));

nFile = numel(AllModelfiles);

for k = 1:nFile
    if strcmp(AllModelfiles(k).name(1:end-10),'FQ_RPE_CK')
        display('1');
        load(fullfile(AllModelfiles(k).folder,AllModelfiles(k).name));
    end

end
fitparMat = cell2mat(fitpar');
% use the median as the animal's "true" parameters
player.params=nanmedian(fitparMat,1);

%% alter different beta 
bkOverSum = 0:0.05:1;
pStayList = zeros(1, length(bkOverSum));
rewardList = zeros(1, length(bkOverSum));
entropyList = zeros(1,length(bkOverSum));
bSum = player.params(2) + player.params(4);

for ii = 1:length(bkOverSum)
    player1.label='FQ_RPE_CK';
    player1.params.a=player.params(1);    % learning rate (also = 1 minus the forgetting rate)
    player1.params.b = bSum*(1-bkOverSum(ii));
    player1.params.ac = player.params(3);
    player1.params.bc = bkOverSum(ii)*bSum;
    %player1.bias = 1;  % if it is a biased payoff matrix
    
    player2.label='algo2';
    player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
    player2.params.trial_history=400;   %trials older than this number are not considered
    
    tlabel=['F-Q-RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
        ' versus algorithm2, n=' num2str(n_stim) ' trials'];
    
    % player_sim.params = [0.5, -1, 1.5, 0.25];
    
    save_path = fullfile(root_path, ['figs_simulation ' player1.label]);
    
    stats_sim=simPennies(player1,player2,n_stim);
    
    % calculate probability of stay, reward rate and entropy
    
    pStay = sum(stats_sim.c(1:end-1,1) == stats_sim.c(2:end,1))/n_stim;
    
    pStayList(ii) = pStay;
    rewardList(ii) = sum(stats_sim.r(:,1))/n_stim;
    
    % calculate entropy
    
    % the possible choice combinations
    choiceBack = 3;
    combos = de2bi([0:2^choiceBack-1],choiceBack);
    combos = 2*(combos - 0.5);  %to make it -1 or 1
    
    % classify animal's choice sequence
    nTrial = size(stats_sim.c,1);
    cumuoccur = zeros(2^choiceBack,nTrial);   %how many times each combo occurred, cumulatively
    for j = choiceBack:nTrial
        c = stats_sim.c(j-choiceBack+1:j,1)';
        idx = ismember(combos,c,'rows');   %find if it matches one of the combos
        if sum(idx)==1
            cumuoccur(:,j) = cumuoccur(:,j-1) + idx;  % update the cumulative occurrence matrix
        else
            cumuoccur(:,j) = cumuoccur(:,j-1);
        end
    end
    p = cumuoccur(:,end)/sum(cumuoccur(:,end));
    entro = -1*sum(p.*log2(p));
    entropyList(ii) = entro;
    %n_plot=500;
    %plot_session_game(stats_sim,n_plot,tlabel);


end

%% simulate the behavior with true value
player1.label='FQ_RPE_CK';
player1.params.a=player.params(1);    % learning rate (also = 1 minus the forgetting rate)
player1.params.b = player.params(2);
player1.params.ac = player.params(3);
player1.params.bc = player.params(4);
%player1.bias = 1;  % if it is a biased payoff matrix

player2.label='algo2';
player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
player2.params.trial_history=400;   %trials older than this number are not considered

tlabel=['F-Q-RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
    ' versus algorithm2, n=' num2str(n_stim) ' trials'];

% player_sim.params = [0.5, -1, 1.5, 0.25];

save_path = fullfile(root_path, ['figs_simulation ' player1.label]);

stats_sim=simPennies(player1,player2,n_stim);

% calculate probability of stay

pStayFit = sum(stats_sim.c(1:end-1,1) == stats_sim.c(2:end,1))/n_stim;
rewardFit = sum(stats_sim.r(:,1))/n_stim;
% the possible choice combinations
choiceBack = 3;
combos = de2bi([0:2^choiceBack-1],choiceBack);
combos = 2*(combos - 0.5);  %to make it -1 or 1

% classify animal's choice sequence
nTrial = size(stats_sim.c,1);
cumuoccur = zeros(2^choiceBack,nTrial);   %how many times each combo occurred, cumulatively
for j = choiceBack:nTrial
    c = stats_sim.c(j-choiceBack+1:j,1)';
    idx = ismember(combos,c,'rows');   %find if it matches one of the combos
    if sum(idx)==1
        cumuoccur(:,j) = cumuoccur(:,j-1) + idx;  % update the cumulative occurrence matrix
    else
        cumuoccur(:,j) = cumuoccur(:,j-1);
    end
end
p = cumuoccur(:,end)/sum(cumuoccur(:,end));
entro = -1*sum(p.*log2(p));
entropyFit= entro;

 %% plot the results
 
fig = figure;
left_color = [0 0 0];
right_color = [1 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
plot(bkOverSum,pStayList,'k');
hold on;
plot(bkOverSum,rewardList,'k--');
hold on;
scatter(player.params(4)/(player.params(2)+player.params(4)), pStayFit,120,'k','filled','HandleVisibility','off');
hold on;
scatter(player.params(4)/(player.params(2)+player.params(4)), rewardFit,120,'k','HandleVisibility','off');
ylim([0 1]);
yyaxis right
plot(bkOverSum, entropyList,'Color',right_color);
ylim([0 3]);
scatter(player.params(4)/(player.params(2)+player.params(4)), entropyFit,120,right_color,'filled','HandleVisibility','off');
xlabel('\beta_K/(\beta+\beta_K)');
lgd = legend('pStay','Reward','Entropy');
set(lgd,'color','none','box','off');

title(['\alpha = ',num2str(player.params(1)),' Sum(\beta) = ',num2str(player.params(2)+player.params(4)),' \alpha_K = ',num2str(player.params(3))]);
set(gca,'box','off');

savepath = fullfile(saveSimPath, 'alter_beta');
print(gcf,'-dpng',savepath );   %png format
saveas(gcf, savepath , 'fig');
saveas(gcf, savepath , 'svg');
%% alter beta_K + beta
% % using median value from fit with animal data
% 
% 
% betaSum = 0.05:0.1:3;
% pStayList = zeros(1, length(betaSum));
% rewardList = zeros(1, length(betaSum));
% entropyList = zeros(1,length(betaSum));
% 
% for ii = 1:length(betaSum)
%     player1.label='FQ_RPE_CK';
%     player1.params.a=player.params(1);    % learning rate (also = 1 minus the forgetting rate)
%     player1.params.b = player.params(2)*betaSum(ii)/(player.params(2)+player.params(4));
%     player1.params.ac = player.params(3);
%     player1.params.bc = player.params(4)*betaSum(ii)/(player.params(2)+player.params(4));
%     %player1.bias = 1;  % if it is a biased payoff matrix
%     
%     player2.label='algo2';
%     player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
%     player2.params.trial_history=400;   %trials older than this number are not considered
%     
%     tlabel=['F-Q-RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
%         ' versus algorithm2, n=' num2str(n_stim) ' trials'];
%     
%     % player_sim.params = [0.5, -1, 1.5, 0.25];
%     
%     save_path = fullfile(root_path, ['figs_simulation ' player1.label]);
%     
%     stats_sim=simPennies(player1,player2,n_stim);
%     
%     % calculate probability of stay, reward rate and entropy
%     
%     pStay = sum(stats_sim.c(1:end-1,1) == stats_sim.c(2:end,1))/n_stim;
%     pStayList(ii) = pStay;
%     rewardList(ii) = sum(stats_sim.r(:,1))/n_stim;
%     
%     % calculate entropy
%     
%     % the possible choice combinations
%     choiceBack = 3;
%     combos = de2bi([0:2^choiceBack-1],choiceBack);
%     combos = 2*(combos - 0.5);  %to make it -1 or 1
%     
%     % classify animal's choice sequence
%     nTrial = size(stats_sim.c,1);
%     cumuoccur = zeros(2^choiceBack,nTrial);   %how many times each combo occurred, cumulatively
%     for j = choiceBack:nTrial
%         c = stats_sim.c(j-choiceBack+1:j,1)';
%         idx = ismember(combos,c,'rows');   %find if it matches one of the combos
%         if sum(idx)==1
%             cumuoccur(:,j) = cumuoccur(:,j-1) + idx;  % update the cumulative occurrence matrix
%         else
%             cumuoccur(:,j) = cumuoccur(:,j-1);
%         end
%     end
%     p = cumuoccur(:,end)/sum(cumuoccur(:,end));
%     entro = -1*sum(p.*log2(p));
%     entropyList(ii) = entro;
%     %n_plot=500;
%     %plot_session_game(stats_sim,n_plot,tlabel);
% 
% 
% end
% figure;
% plot(betaSum,pStayList,'k');
% hold on;
% plot(betaSum,rewardList,'r');
% hold on;
% plot(betaSum, entropyList,'b');
% hold on;
% scatter(player.params(2)+player.params(4), pStayFit,120,'k','filled');
% hold on;
% scatter(player.params(2)+player.params(4), rewardFit,120,'r','filled');
% hold on;
% scatter(player.params(2)+player.params(4), entropyFit,120,'b','filled');
% xlabel('\beta+\beta_K');
% legend('pStay','reward');
% title(['\alpha = ',num2str(player.params(1)),' \beta = ',num2str(player.params(2)),' \alpha_K = ',num2str(player.params(3))]);
% 
% savepath = fullfile(saveSimPath, 'alter_betaSum');
% print(gcf,'-dpng',savepath );   %png format
% saveas(gcf, savepath , 'fig');
% 
% % what if the betaK/(betaK+beta) is less than 0.5?
% betaSum = 0.05:0.1:3;
% pStayList = zeros(1, length(betaSum));
% rewardList = zeros(1, length(betaSum));
% entropyList = zeros(1,length(betaSum));
% 
% for ii = 1:length(betaSum)
%     player1.label='FQ_RPE_CK';
%     player1.params.a=player.params(1);    % learning rate (also = 1 minus the forgetting rate)
%     player1.params.b = 0.25*betaSum(ii);
%     player1.params.ac = player.params(3);
%     player1.params.bc = 0.25*betaSum(ii);
%     %player1.bias = 1;  % if it is a biased payoff matrix
%     
%     player2.label='algo2';
%     player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
%     player2.params.trial_history=400;   %trials older than this number are not considered
%     
%     tlabel=['F-Q-RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
%         ' versus algorithm2, n=' num2str(n_stim) ' trials'];
%     
%     % player_sim.params = [0.5, -1, 1.5, 0.25];
%     
%     save_path = fullfile(root_path, ['figs_simulation ' player1.label]);
%     
%     stats_sim=simPennies(player1,player2,n_stim);
%     
%     % calculate probability of stay, reward rate and entropy
%     
%     pStay = sum(stats_sim.c(1:end-1,1) == stats_sim.c(2:end,1))/n_stim;
%     pStayList(ii) = pStay;
%     rewardList(ii) = sum(stats_sim.r(:,1))/n_stim;
%     
%     % calculate entropy
%     
%     % the possible choice combinations
%     choiceBack = 3;
%     combos = de2bi([0:2^choiceBack-1],choiceBack);
%     combos = 2*(combos - 0.5);  %to make it -1 or 1
%     
%     % classify animal's choice sequence
%     nTrial = size(stats_sim.c,1);
%     cumuoccur = zeros(2^choiceBack,nTrial);   %how many times each combo occurred, cumulatively
%     for j = choiceBack:nTrial
%         c = stats_sim.c(j-choiceBack+1:j,1)';
%         idx = ismember(combos,c,'rows');   %find if it matches one of the combos
%         if sum(idx)==1
%             cumuoccur(:,j) = cumuoccur(:,j-1) + idx;  % update the cumulative occurrence matrix
%         else
%             cumuoccur(:,j) = cumuoccur(:,j-1);
%         end
%     end
%     p = cumuoccur(:,end)/sum(cumuoccur(:,end));
%     entro = -1*sum(p.*log2(p));
%     entropyList(ii) = entro;
%     %n_plot=500;
%     %plot_session_game(stats_sim,n_plot,tlabel);
% 
% 
% end
% figure;
% plot(betaSum,pStayList,'k');
% hold on;
% plot(betaSum,rewardList,'r');
% hold on;
% plot(betaSum, entropyList,'b');
% hold on;
% scatter(player.params(2)+player.params(4), pStayFit,120,'k','filled');
% hold on;
% scatter(player.params(2)+player.params(4), rewardFit,120,'r','filled');
% hold on;
% scatter(player.params(2)+player.params(4), entropyFit,120,'b','filled');
% xlabel('\beta+\beta_K');
% legend('pStay','reward');
% title(['\alpha = ',num2str(player.params(1)),' \beta = ',num2str(player.params(2)),' \alpha_K = ',num2str(player.params(3))]);
% 
% savepath = fullfile(saveSimPath, 'alter_betaSum_betaK0.25');
% print(gcf,'-dpng',savepath );   %png format
% saveas(gcf, savepath , 'fig');
% 
% %% now alter alpha
% akOvera = 0.1:0.1:2;
% pStayList = zeros(1, length(akOvera));
% rewardList = zeros(1, length(akOvera));
% entropyList = zeros(1, length(akOvera));
% 
% for ii = 1:length(akOvera)
%     player1.label='FQ_RPE_CK';
%     player1.params.a=player.params(1);    % learning rate (also = 1 minus the forgetting rate)
%     player1.params.b = player.params(2);
%     player1.params.ac = player.params(1)*akOvera(ii);
%     player1.params.bc = player.params(4);
%     %player1.bias = 1;  % if it is a biased payoff matrix
%     
%     player2.label='algo2';
%     player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
%     player2.params.trial_history=400;   %trials older than this number are not considered
%     
%     tlabel=['F-Q-RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
%         ' versus algorithm2, n=' num2str(n_stim) ' trials'];
%     
%     % player_sim.params = [0.5, -1, 1.5, 0.25];
%     
%     save_path = fullfile(root_path, ['figs_simulation ' player1.label]);
%     
%     stats_sim=simPennies(player1,player2,n_stim);
%     
%     % calculate probability of stay
%     
%     pStay = sum(stats_sim.c(1:end-1,1) == stats_sim.c(2:end,1))/n_stim;
%     pStayList(ii) = pStay;
%     rewardList(ii) = sum(stats_sim.r(:,1))/n_stim;
%     n_plot=500;
%     %plot_session_game(stats_sim,n_plot,tlabel);
%     % calculate entropy
%     
%     % the possible choice combinations
%     choiceBack = 3;
%     combos = de2bi([0:2^choiceBack-1],choiceBack);
%     combos = 2*(combos - 0.5);  %to make it -1 or 1
%     
%     % classify animal's choice sequence
%     nTrial = size(stats_sim.c,1);
%     cumuoccur = zeros(2^choiceBack,nTrial);   %how many times each combo occurred, cumulatively
%     for j = choiceBack:nTrial
%         c = stats_sim.c(j-choiceBack+1:j,1)';
%         idx = ismember(combos,c,'rows');   %find if it matches one of the combos
%         if sum(idx)==1
%             cumuoccur(:,j) = cumuoccur(:,j-1) + idx;  % update the cumulative occurrence matrix
%         else
%             cumuoccur(:,j) = cumuoccur(:,j-1);
%         end
%     end
%     p = cumuoccur(:,end)/sum(cumuoccur(:,end));
%     entro = -1*sum(p.*log2(p));
%     entropyList(ii) = entro;
%     %n_plot=500;
% 
% 
% end
%  %% plot the results
% 
% 
% 
% 
% 
% fig = figure;
% left_color = [0 0 0];
% right_color = [1 0 1];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% 
% yyaxis left
% 
% plot(akOvera,pStayList,'k');
% hold on;
% plot(akOvera,rewardList,'k--');
% hold on;
% 
% scatter(player.params(3)/(player.params(1)), pStayFit,120,'k','filled','HandleVisibility','off');
% hold on;
% scatter(player.params(3)/(player.params(1)), rewardFit,120,'k','HandleVisibility','off');
% ylim([0 1]);
% 
% yyaxis right
% 
% plot(akOvera, entropyList,'b');
% hold on;
% scatter(player.params(3)/(player.params(1)), entropyFit,120,right_color,'filled','HandleVisibility','off');
% ylim([0 3]);
% xlabel('\alpha_K/\alpha');
% lgd = legend('pStay','Reward','Entropy');
% set(lgd,'color','none','box','off');
% title(['\alpha = ',num2str(player.params(1)),' \beta = ',num2str(player.params(2)),' \beta_K = ',num2str(player.params(4))]);
% 
% savepath = fullfile(saveSimPath, 'alter_alpha');
% print(gcf,'-dpng',savepath );   %png format
% saveas(gcf, savepath , 'fig');
% saveas(gcf, savepath , 'svg');
% 
% 

