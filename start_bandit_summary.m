%% Summarize multiple sessions of behavioral performance
% run this after running start_bandit

clearvars;
close all;

%% setup path and plotting formats

value_setPathList;

setup_figprop;  %set up default figure plotting parameters

%% load data file list

m=2;

switch m
    case 1  % Load the data set from the pilot bandit data
        data_subdir = fullfile(data_dir,'pilot bandit');
        [ dirs, expData ] = expData_bandit(data_subdir);
    case 2  % Load the data set from Wanyu's bandit data, fixed inter-trial interval
        data_subdir = fullfile(data_dir,'pilot bandit-fixedITI');
        [ dirs, expData ] = expData_banditFixed(data_subdir);
end

%% load the data

stats_all.c=[];  %concatenate choices and outcomes across sessions
stats_all.r=[];
for i = 1:numel(expData)

    if isfield(expData(i),'onefolder')   %all the log files saved in one folder
        temp = sprintf('%04d',expData(i).logfilenum);
        cd(fullfile(dirs.analysis,expData(i).sub_dir,temp));
    else
        cd(fullfile(dirs.analysis,expData(i).sub_dir));
    end
    load('beh.mat');

    lick_trType_array{i}=lick_trType;
    
    iti_array{i}=iti_trType;
    respTime_array{i}=respTime_trType;
    
    sw_array{i}=sw_output;
    sw_hrside_array{i}=sw_hrside_output;
    bl_array{i}=bl_output;
    
    lregRCUC_array{i}=lregRCUC_output;
    lregCRInt_array{i}=lregCRInt_output;

    if exist('nlike_array','var')
        fname = fieldnames(nlike);
        for j=1:numel(fname)  %append for each field
            nlike_array.(fname{j})=[nlike_array.(fname{j}); nlike.(fname{j})];
            bic_array.(fname{j})=[bic_array.(fname{j}); bic.(fname{j})];
        end
    else
        nlike_array = nlike;
        bic_array = bic;
    end
    
    stats_all.c=[stats_all.c; stats.c];
    stats_all.r=[stats_all.r; stats.r];
    
    close all;
    clearvars -except i dirs expData ...
        lick_trType_array lregRCUC_array lregCRInt_array...
        sw_array sw_hrside_array bl_array...
        iti_array respTime_array nlike_array bic_array...
        stats_all;
end

%%
savebehfigpath = fullfile(dirs.summary,'figs-beh');
if ~exist(savebehfigpath,'dir')
    mkdir(savebehfigpath);
end

cd(savebehfigpath);
tlabel=strcat('Group summary, n=',int2str(numel(iti_array)));

%%
plot_lickrate_byTrialType(lick_trType_array);

plot_val_byTrialType(respTime_array);
print(gcf,'-dpng','rt_byTrialType');    %png format
saveas(gcf, 'rt_byTrialType', 'fig');

plot_val_byTrialType(iti_array);
print(gcf,'-dpng','iti_byTrialType');   %png format
saveas(gcf, 'iti_byTrialType', 'fig');

%% plot choice behavior - around switches and during block

rule_labels = {'0.7:0.1','0.1:0.7'};
plot_switch(sw_array,tlabel,rule_labels);

plot_switch_hrside(sw_hrside_array,tlabel);

plot_block(bl_array,tlabel,rule_labels);
    
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
x=1;    %analyze player 1

%regressors = rewarded choice, unrewarded choice
[lregRCUC_output, ~, ~, ~]=logreg_RCUC(stats_all,1,num_regressor);
plot_logreg(lregRCUC_output,tlabel);
print(gcf,'-dpng','logregRCUC_concat');    %png format
saveas(gcf, 'logregRCUC_concat', 'fig');

%regressors = choice, choice x reward
[lregCRInt_output, ~, ~, ~]=logreg_CRInt(stats_all,1,num_regressor);
plot_logreg(lregCRInt_output,tlabel);
print(gcf,'-dpng','logregCRInt_concat');    %png format
saveas(gcf, 'logregCRInt_concat', 'fig');

%% fit to models
disp(['Total number of trials: ' int2str(sum(~isnan(stats_all.c(:,1))))]);
disp(['Model fit, using concatenated choice behaviors:']);

fun = 'WSLSfun';
initpar=0.5971; % initial [prob_WSLS]
lb=0;
ub=1;
[fitpar.WSLS, ~, bic.WSLS, nlike.WSLS]=fit_fun(stats_all,fun,initpar,1,lb,ub);

fun = 'Q_RPEfun';
initpar=[0.67 2.26]; % initial [alpha beta]
[fitpar.Q_RPE, ~, bic.Q_RPE, nlike.Q_RPE]=fit_fun(stats_all,fun,initpar,1);

fun = 'DFQfun';
initpar=[0.70 0.30 3.51 -1.53]; % initial [alpha1 alpha2 kappa1 kappa2]
[fitpar.DFQ, ~, bic.DFQ, nlike.DFQ]=fit_fun(stats_all,fun,initpar,1);

fun = 'FQfun';
initpar=[0.52 3.68 -1.13]; % initial [alpha1 kappa1 kappa2]
[fitpar.FQ, ~, bic.FQ, nlike.FQ]=fit_fun(stats_all,fun,initpar,1);

fun = 'Qfun';
initpar=[0.66 5.00 -2.74]; % initial [alpha1 kappa1 kappa2]
[fitpar.Q, ~, bic.Q, nlike.Q]=fit_fun(stats_all,fun,initpar,1);

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

%% plays sound when done
load train;
sound(y,Fs);
