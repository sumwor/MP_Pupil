%% Summarize multiple sessions of behavioral performance
% run this after running start_pennies

clearvars;
close all;

%% setup path and plotting formats

value_setPathList;

setup_figprop;  %set up default figure plotting parameters

%% load data file list

m=1;

switch m
    case 1  % Load the data set from the pilot MP data
        data_subdir = fullfile(data_dir);
        [ dirs, expData ] = expData_pennies(data_subdir);
end

%% load the data

stats_all.c=[];  %concatenate choices and outcomes across sessions
stats_all.r=[];
% subMask = [];
for i = 1:numel(expData)
%     switch expData(i).logfile(1:3)
%         case '862'
%             mask = 1;
%         case '863'
%             mask = 2;
%         case '866'
%             mask = 3;
%         case '867'
%             mask = 4;
%         case '868'
%             mask = 5;
%         case '880'
%             mask = 6;
%         case '882'
%             mask = 7;
%         case '883'
%             mask = 8;
%     end
%     subMask = [subMask, mask];
    if isfield(expData(i),'onefolder')   %all the log files saved in one folder
        temp = sprintf('%04d',expData(i).logfilenum);
        cd(fullfile(dirs.analysis,expData(i).sub_dir,temp));
    else
        cd(fullfile(dirs.analysis,expData(i).sub_dir));
    end
    load('beh.mat');

    lick_trType_array{i}=lick_trType;
    
    lregRCUC_array{i}=lregRCUC_output;
    lregCRInt_array{i}=lregCRInt_output;
    
    iti_array{i}=iti_trType;
    respTime_array{i}=respTime_trType;
    trueRespTime{i} = trialData.rt;
    choiceBySession{i} = stats;
    
    if exist('nlike_array','var')
        fname = fieldnames(nlike);
        for j=1:numel(fname)  %append for each field
            nlike_array.(fname{j})=[nlike_array.(fname{j}); nlike.(fname{j})];
            bic_array.(fname{j})=[bic_array.(fname{j}); bic.(fname{j})];
        end
        fname = fieldnames(fitpar);
        for j=1:numel(fname)  %append for each field
            fitpar_array.(fname{j})=[fitpar_array.(fname{j}); fitpar.(fname{j})];
        end
    else
        nlike_array = nlike;
        bic_array = bic;
        fitpar_array = fitpar;
    end
    
    nTrial_array(i)=sum(stats.c(:,1)==-1)+sum(stats.c(:,1)==1);
    entro_array(i)=entro;
    rrate_array(i)=sum(stats.r==1)/(sum(stats.r==1)+sum(stats.r==0));
    
    stats_all.c=[stats_all.c; stats.c];
    stats_all.r=[stats_all.r; stats.r];
    
    close all;
    clearvars -except i dirs expData ...
        lick_trType_array lregRCUC_array lregCRInt_array iti_array respTime_array trueRespTime choiceBySession...
        nlike_array bic_array fitpar_array...
        nTrial_array entro_array rrate_array subMask...
        stats_all;
end


%%
savebehfigpath = fullfile(dirs.summary,'figs-beh');
if ~exist(savebehfigpath,'dir')
    mkdir(savebehfigpath);
end

cd(savebehfigpath);
%save('stats_1.mat', 'choiceBySession','entro_array','respTime_array', 'rrate_array', 'subMask');
tlabel=strcat('Group summary, n=',int2str(numel(iti_array)));

%%
plot_lickrate_byTrialType(lick_trType_array);

plot_val_byTrialType(respTime_array);
print(gcf,'-dpng','rt_byTrialType');    %png format
saveas(gcf, 'rt_byTrialType', 'fig');

plot_val_byTrialType(iti_array);
print(gcf,'-dpng','iti_byTrialType');   %png format
saveas(gcf, 'iti_byTrialType', 'fig');

%% fit with logistic regression

%average of sessions
plot_logreg(lregRCUC_array,tlabel);
print(gcf,'-dpng','logregRCUC');    %png format
saveas(gcf, 'logregRCUC', 'fig');

plot_logreg(lregCRInt_array,tlabel);
print(gcf,'-dpng','logregCRInt');    %png format
saveas(gcf, 'logregCRInt', 'fig');

%% concatenate sessions for each animal, fit with logistic regression 
num_regressor=15;

%regressors = rewarded choice, unrewarded choice
[lregRCUC_output, ~, ~, ~]=logreg_RCUC(stats_all,1,num_regressor);
plot_logreg(lregRCUC_output,tlabel);
print(gcf,'-dpng','logregRCUC_concat');    %png format
saveas(gcf, 'logregRCUC_concat', 'fig');

%regressors = choice, choice x reward (equivalent to computer's choice)
[lregCRInt_output, ~, ~, ~]=logreg_CRInt(stats_all,1,num_regressor);
plot_logreg(lregCRInt_output,tlabel);
print(gcf,'-dpng','logregCRInt_concat');    %png format
saveas(gcf, 'logregCRInt_concat', 'fig');

%%
disp(['Total number of trials: ' int2str(sum(~isnan(stats_all.c(:,1))))]);
disp(['Model fit, using concatenated choice behaviors:']);

fun = 'WSLSfun';
initpar=0.5344; % initial [prob_WSLS]
lb=0;
ub=1;
[fitpar.WSLS, ~, bic.WSLS, nlike.WSLS]=fit_fun(stats_all,fun,initpar,1,lb,ub);

fun = 'Q_RPEfun';
initpar=[1.1507 0.1759]; % initial [alpha beta]
[fitpar.Q_RPE, ~, bic.Q_RPE, nlike.Q_RPE]=fit_fun(stats_all,fun,initpar,1);

fun = 'DFQfun';
initpar=[0.118 0.0147 3.7881 -3.1804]; % initial [alpha1 alpha2 kappa1 kappa2]
[fitpar.DFQ, ~, bic.DFQ, nlike.DFQ]=fit_fun(stats_all,fun,initpar,1);

fun = 'FQfun';
initpar=[0.06 2.01 -1.04]; % initial [alpha1 kappa1 kappa2]
[fitpar.FQ, ~, bic.FQ, nlike.FQ]=fit_fun(stats_all,fun,initpar,1);

fun = 'Qfun';
initpar=[1.13 -2.12 2.29]; % initial [alpha1 kappa1 kappa2]
[fitpar.Q, ~, bic.Q, nlike.Q]=fit_fun(stats_all,fun,initpar,1);

fun = 'FQfun_withbeta';
initpar=[0.1 0 1]; % initial [alpha1 kappa2 beta]
[fitpar.FQ_wB, ~, bic.FQ_wB, nlike.FQ_wB]=fit_fun(stats_all,fun,initpar,1);

fun = 'FQfun_withbetabias';
initpar=[0.1 0 1 0]; % initial [alpha1 kappa2 beta deltaQ]
[fitpar.FQ_wBB, ~, bic.FQ_wBB, nlike.FQ_wBB]=fit_fun(stats_all,fun,initpar,1);

fun = 'FQfun_withbetaCA';
initpar=[0.1 0 1 0 0]; % initial [alpha1 kappa2 beta tau phi]
[fitpar.FQ_wBC, ~, bic.FQ_wBC, nlike.FQ_wBC]=fit_fun(stats_all,fun,initpar,1);

[~, ~, bic.logregCRInt2, nlike.logregCRInt2]=logreg_CRInt(stats_all,1,2);
[~, ~, bic.logregCRInt5, nlike.logregCRInt5]=logreg_CRInt(stats_all,1,5);
[~, ~, bic.logregCRInt10, nlike.logregCRInt10]=logreg_CRInt(stats_all,1,10);

fitpar
bic

%% simulate to get latent action value, compare with animal's choice
player1.label='algo_DFQ';
player1.params.a1=fitpar.FQ(1);
player1.params.a2=fitpar.FQ(1);
player1.params.k1=fitpar.FQ(2);
player1.params.k2=fitpar.FQ(3);

stats_sim=predictAgent(player1,stats_all);

x = 1;  %player 1
n_plot = 500;   %plot first 500 trials
plot_session_qparam(stats_sim,x,n_plot);

%%
disp('---Mean normalized likelihood for fits per session');
fname = fieldnames(nlike_array);
for j=1:numel(fname)  %append for each field
    disp([fname{j} ' - ' num2str(nanmean(nlike_array.(fname{j})))]);
end

disp('---Mean BIC for fits per session');
fname = fieldnames(bic_array);
for j=1:numel(fname)  %append for each field
    disp([fname{j} ' - ' num2str(nanmean(bic_array.(fname{j})))]);
end

%%
figure;
subplot(2,5,1); hold on;
plot(rand(1,numel(nTrial_array)),nTrial_array,'k^','MarkerSize',15);
ylim([0 550]); xlim([-1 2]);
ylabel('Trials performed');

subplot(2,5,3); hold on;
plot(rand(1,numel(entro_array)),entro_array,'k^','MarkerSize',15);
plot([-1 2],[3 3],'k--','LineWidth',2);
ylim([2.5 3.1]); xlim([-1 2]);
ylabel('Entropy (bits)');

subplot(2,5,5); hold on;
plot(rand(1,numel(rrate_array)),100*rrate_array,'k^','MarkerSize',15);
plot([-1 2],[50 50],'k--','LineWidth',2);
ylim([35 65]); xlim([-1 2]);
ylabel('Reward rate (%)');

subplot(2,3,4); hold on;
plot([-0.1 1.1],[0 0],'k--','LineWidth',2);
plot([0 0],[-0.4 4],'k--','LineWidth',2);
plot(fitpar_array.Q_RPE(:,1),fitpar_array.Q_RPE(:,2),'k.','MarkerSize',30);
xlabel('\alpha, Learning rate');
ylabel('\beta, Inverse temperature');
xlim([-0.1 1.1]);
ylim([-0.4 4]);
print(gcf,'-dpng',['alpha-beta' int2str(x)]);    %png format
saveas(gcf,['alpha-beta' int2str(x)], 'fig');

alpha = fitpar_array.Q_RPE(:,1); beta = fitpar_array.Q_RPE(:,2);
save('a_b_saline.mat','alpha','beta');

%% response time - total value correlation 

% use the estimation of DFQ algorithm to get the estimation of Ql and Qr

% numAnimal = expData(end).animalID;
qSum = cell(0);
for jj = 1:length(expData)
    player1.label='algo_DFQ';
    player1.params.a1=fitpar_array.DFQ(jj, 1);
    player1.params.a2=fitpar_array.DFQ(jj, 2);
    player1.params.k1=fitpar_array.DFQ(jj, 3);
    player1.params.k2=fitpar_array.DFQ(jj, 4);
    
    stats_sim=predictAgent(player1,choiceBySession{jj});
    qSum{jj} = stats_sim.ql + stats_sim.qr;
    qDif{jj} = stats_sim.ql - stats_sim.qr;
    if mod(jj, 10) == 0
        jj
    end

end

% correlate with response time
% get rid of the NaN
maxLag = zeros(1, length(qSum));
maxCorr = zeros(1, length(qSum));
corrZero = zeros(1, length(qSum));
corrSum = zeros(1, length(qSum));
pCorrSum = zeros(1, length(qSum));
for tt = 1:length(qSum)
    noNan=~isnan(trueRespTime{tt});
    [R,P] = corrcoef(abs(qSum{tt}(noNan)), trueRespTime{tt}(noNan)); 
    corrSum(tt) = R(2,1);
    pCorrSum(tt) = P(2,1);
%     noNan=~isnan(trueRespTime{tt});
%     [c,lags] = xcorr(qSum{tt}(noNan), trueRespTime{tt}(noNan), 'coeff'); 
%     [maxcorr, ind] = max(abs(c));
%     Lags{tt} = lags;
%     Corr{tt} = c;
%     maxLag(tt) = lags(ind);
%     maxCorr(tt) = c(ind);
%     corrZero(tt) = xcorr(qSum{tt}(noNan), trueRespTime{tt}(noNan),0, 'coeff'); 
% %     figure;
%     stem(lags,c)
end
figure;
subplot(2,1,1);histogram(maxCorr);

binSize = 0.1;
corrBin_total = zeros(1, length([-1:0.1:1]));
corrBin_sigtotal = zeros(1, length([-1:0.1:1]));
for ii = 1:length(corrSum)
    corrBin_total(ceil((corrSum(ii)+1)/0.1)) = corrBin_total(ceil((corrSum(ii)+1)/0.1)) + 1;
    if pCorrSum(ii) < 0.05
        corrBin_sigtotal(ceil((corrSum(ii)+1)/0.1)) =  corrBin_sigtotal(ceil((corrSum(ii)+1)/0.1)) + 1;
    end
end
figure;
bar([-1:0.1:1], corrBin_total, 'FaceColor', 'black', 'FaceAlpha',0.3);
hold on; bar(-1:0.1:1, corrBin_sigtotal, 'FaceColor','Black','FaceAlpha',0.8);
title('Correlations between total value and response time');
xlabel('Correlation coefficients');
ylabel('# of sessions')
print(gcf,'-dpng','pos_neg_sessions_total');    %png format
saveas(gcf, 'pos_neg_sessions_total', 'fig');

% figure;
% histogram(corrZero, 10, 'FaceColor', 'black');
% title('Histogram correlations between total value and response time');
% xlabel('Subjects');
% ylabel('# of sessions')
% print(gcf,'-dpng','pos_neg_sessions');    %png format
% saveas(gcf, 'pos_neg_sessions', 'fig');
% get the below -.5 or above .5 sessions
wholeInd = 1:length(corrZero);
% sessionPos = corrS > 0.5;
% sessionNeg = corrZero < -0.5;
sessionPos = corrSum > 0;
sessionNeg = corrZero < 0;
sessionPos_sig = pCorrSum < 0.05 & corrSum > 0;
sessionNeg_sig = pCorrSum < 0.05 & corrSum <0;

% get identity mask
IdMask = zeros(1, length(corrZero));
for ll = 1:length(corrZero)
    IdMask(ll) = expData(ll).animalID;
end
% get the session identity
% make bin counts myself
pos_total = zeros(1,12); neg_total = zeros(1,12); 
pos_sigtotal = zeros(1,12); neg_sigtotal = zeros(1,12);
for ii = 1:12
    pos_total(ii) = sum(IdMask(sessionPos)==ii);
    neg_total(ii) = sum(IdMask(sessionNeg)==ii);
    pos_sigtotal(ii) = sum(IdMask(sessionPos_sig)==ii);
    neg_sigtotal(ii) = sum(IdMask(sessionNeg_sig)==ii);
end

figure; bar([pos_total;neg_total]', 1,'FaceColor','flat','EdgeColor', 'flat' ,'EdgeAlpha',0.5, 'FaceAlpha', 0.5);
%hold on; bar([1.5:12.5],neg, 'FaceColor', [0.8500 0.3250 0.0980],'FaceAlpha', 0.3);
%hold on; bar([1:12], pos_sig, 'FaceColor', [0.8500 0.3250 0.0980],'FaceAlpha', 0.8);
hold on; bar([pos_sigtotal; neg_sigtotal]', 1,'FaceColor','flat' ,'EdgeAlpha',0.8,  'FaceAlpha', 0.8);
legend('pos correlation', 'neg correlation');
title('Positive and negative correlation between total value and response time by subjects');
xlabel('Subjects');
ylabel('# of sessions')
print(gcf,'-dpng','pos_neg_sessions_subjects_total');    %png format
saveas(gcf, 'pos_neg_sessions_subjects_total', 'fig');

% figure; histogram(IdMask(sessionPos),12);
% hold on; histogram(IdMask(sessionNeg), 12);
% legend('pos correlation', 'neg correlation');
% title('Positive and negative correlation between total value and response time by subjects');
% xlabel('Subjects');
% ylabel('# of sessions')
% print(gcf,'-dpng','pos_neg_sessions_subjects');    %png format
% saveas(gcf, 'pos_neg_sessions_subjects', 'fig');

% test the total value
[upLimit, upInd] = max(maxLag);
[lowerLimit, lowInd] = min(maxLag);
% plot the highest and lowest data
figure;
subplot(2,1,1);
stem(Lags{upInd}, Corr{upInd});
subplot(2,1,2);
stem(Lags{lowInd}, Corr{lowInd});

interval = floor(lowerLimit/10)*10 - 5:10:ceil(upLimit/10)*10 + 5;
lagCount = zeros(1, length(interval) -1 );
for kk = 1:length(maxLag)
    lagCount(floor((maxLag(kk) - interval(1)) / 10) + 1) = lagCount(floor((maxLag(kk) - interval(1)) / 10) + 1) + 1;
end
figure;plot(interval(1:end-1),lagCount/sum(lagCount));

% normalize (zscore)
for zz = 1:length(trueRespTime)
    normtrueRespTime{zz} = (trueRespTime{zz} - nanmean(trueRespTime{zz})) / nanstd(trueRespTime{zz});
    normQSum{zz} = (qSum{zz} - mean(qSum{zz})) / std(qSum{zz});
    normQDiff{zz} = (qDif{zz} - mean(qDif{zz})) / std(qDif{zz});
end

% test the value difference and the response time

corrDif = zeros(1, length(qDif));
pCorrDif = zeros(1, length(qDif));
for tt = 1:length(qDif)
    noNan=~isnan(trueRespTime{tt});
    [R,P] = corrcoef(abs(qDif{tt}(noNan)), trueRespTime{tt}(noNan)); 
    corrDif(tt) = R(2,1);
    pCorrDif(tt) = P(2,1);
end

sessionPos_dif = corrDif > 0;
sessionNeg_dif = corrDif < -0;
sessionPos_difsig = corrDif > 0 & pCorrDif < 0.05;
sessionNeg_difsig = corrDif < 0 & pCorrDif < 0.05;

% make bin counts myself
pos = zeros(1,12); neg = zeros(1,12); 
pos_sig = zeros(1,12); neg_sig = zeros(1,12);
for ii = 1:12
    pos(ii) = sum(IdMask(sessionPos_dif)==ii);
    neg(ii) = sum(IdMask(sessionNeg_dif)==ii);
    pos_sig(ii) = sum(IdMask(sessionPos_difsig)==ii);
    neg_sig(ii) = sum(IdMask(sessionNeg_difsig)==ii);
end
figure; bar([pos;neg]', 1,'FaceColor','flat','EdgeColor', 'flat' ,'EdgeAlpha',0.5, 'FaceAlpha', 0.5);
%hold on; bar([1.5:12.5],neg, 'FaceColor', [0.8500 0.3250 0.0980],'FaceAlpha', 0.3);
%hold on; bar([1:12], pos_sig, 'FaceColor', [0.8500 0.3250 0.0980],'FaceAlpha', 0.8);
hold on; bar([pos_sig; neg_sig]', 1,'FaceColor','flat' ,'EdgeAlpha',0.8,  'FaceAlpha', 0.8);
legend('pos correlation', 'neg correlation');
title('Positive and negative correlation between absolute value difference and response time by subjects');
xlabel('Subjects');
ylabel('# of sessions')
print(gcf,'-dpng','pos_neg_sessions_subjects_dif');    %png format
saveas(gcf, 'pos_neg_sessions_subjects_dif', 'fig');


binSize = 0.1;
corrBin = zeros(1, length([-1:0.1:1]));
corrBin_sig = zeros(1, length([-1:0.1:1]));
for ii = 1:length(corrDif)
    corrBin(ceil((corrDif(ii)+1)/0.1)) = corrBin(ceil((corrDif(ii)+1)/0.1)) + 1;
    if pCorrDif(ii) < 0.05
        corrBin_sig(ceil((corrDif(ii)+1)/0.1)) =  corrBin_sig(ceil((corrDif(ii)+1)/0.1)) + 1;
    end
end
figure;
bar([-1:0.1:1], corrBin, 'FaceColor', 'black', 'FaceAlpha',0.3);
hold on; bar(-1:0.1:1, corrBin_sig, 'FaceColor','Black','FaceAlpha',0.8);
title('Correlations between absolute value difference and response time');
xlabel('Correlation coefficients');
ylabel('# of sessions')
print(gcf,'-dpng','pos_neg_sessions_dif');    %png format
saveas(gcf, 'pos_neg_sessions_dif', 'fig');


%% check response time and session stage
% first 50 trials ans last 50 trials
sessionLength = cellfun(@length,trueRespTime);
meanRTStart = zeros(1, sum(sessionLength>300));
meanRTEnd = zeros(1, sum(sessionLength>300));
oInd = 1:length(trueRespTime);
Ind = oInd(sessionLength>300);
for ii=1:length(Ind)
    rt = trueRespTime{Ind(ii)}(~isnan(trueRespTime{Ind(ii)}));
    meanRTStart(ii) = mean(rt(1:50));
    meanRTEnd(ii) = mean(rt(end-49:end));
end

[h1,p1] = ttest(meanRTStart, meanRTEnd, 'tail','left')
meanStart = mean(meanRTStart);
meanEnd = mean(meanRTEnd);
stdStart =std(meanRTStart);
stdEnd = std(meanRTEnd);

diff = meanRTEnd - meanRTStart;
figure;boxplot(diff, 'Notch','on');
xlabel('End of session - start of session');
ylabel('Response time difference (s)');
title('Response time difference between session end and session start');
print(gcf,'-dpng','rt_box');    %png format
saveas(gcf, 'rt_box', 'fig');

figure;
bar([meanStart, meanEnd],0.6,'black');
set(gca,'xticklabel', {'First 50 trials', 'Last 50 trials'});
hold on;
er = errorbar([meanStart,meanEnd], [stdStart,stdEnd]);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylim([0 0.5])
title('Averaged response time');
print(gcf,'-dpng','rt_stage');    %png format
saveas(gcf, 'st_stage', 'fig');
%% plays sound when done
load train;
sound(y,Fs);
