%% Analyze performance of a matching pennies game

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

%% process data files
for i = 1:numel(expData)
    disp(['Processing file ' int2str(i) ' out of ' int2str(numel(expData)) '.']);

    % setup/create subdirectories to save analysis and figures
    if isfield(expData(i),'onefolder')   %all the log files saved in one folder
        temp = sprintf('%04d',expData(i).logfilenum);
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir,temp);
        savebehfigpath = fullfile(dirs.analysis,expData(i).sub_dir,temp,'figs-beh');
    else
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
        savebehfigpath = fullfile(dirs.analysis,expData(i).sub_dir,'figs-beh');
    end

    if ~exist(savematpath,'dir')
        mkdir(savematpath);
    end
    if ~exist(savebehfigpath,'dir')
        mkdir(savebehfigpath);
    end

    [ logData ] = MP_parseLogfile( fullfile(dirs.data,expData(i).sub_dir), expData(i).logfile );
    
    [ sessionData, trialData ] = MP_getSessionData( logData );
    
    [ trials ] = MP_getTrialMasks( trialData );
    
    %% setting up data for plotting
    
    stats.playerlabel{1} = 'Mouse';
    
    %reward probabilities for left and right sides
    
    stats.rule=nan(sessionData.nTrials,1);
    stats.rule(trials.L70R10,1)=1;
    stats.rule(trials.L10R70,1)=2;
    rule_labels = {'0.7:0.1','0.1:0.7'};
    probList=[0.7 0.1; 0.1 0.7];
    
    stats.rewardprob=nan(sessionData.nTrials,2);
    for j = 1:size(probList,1)      %number of reward probabilities sets
        for k = 1:size(probList,2)  %left and right
            stats.rewardprob(stats.rule==j,k) = probList(j,k);
        end
    end
    
    %choice: left=-1; right=1; miss=NaN
    stats.c=nan(sessionData.nTrials,1);
    stats.c(trials.left,1)=-1;
    stats.c(trials.right,1)=1;
 
    %outcome: reward=1; no reward:0; miss=NaN
    stats.r=nan(sessionData.nTrials,1);
    stats.r(trials.reward)=1;
    stats.r(trials.noreward)=0;
    
    %% analysis of behavioral performance
    
    % plot behavior in raw format
    tlabel=strcat('Subject=',char(logData.subject),', Time=',char(logData.dateTime(1)),'-',char(logData.dateTime(2)));
    time_range=[-2 6];
%    plot_session_beh_vert(trialData,trials,[],tlabel,time_range);

    % plot choice behavior - whole sessions
    cd(savebehfigpath);
    tlabel=strcat('Subject=',logData.subject,', Time=',logData.dateTime(1),'-',logData.dateTime(2));
    
    plot_session_task(stats,sessionData.nTrials,tlabel);
    
    %% plot choice behavior - around switches
    
    trials_back=10;  % set number of previous trials
    
    sw_output=choice_switch(stats,trials_back);
    plot_switch(sw_output,tlabel,rule_labels);
    
    sw_hrside_output=choice_switch_hrside(stats,trials_back);
    plot_switch_hrside(sw_hrside_output,tlabel);
    
    %% plot choice behavior - during a block
    
    trials_forw=10;  % set number of trials
    
    bl_output=choice_block(stats,trials_forw);
    plot_block(bl_output,tlabel,rule_labels);

    %% plot lick rates
    trialType={'reward','noreward'};
    edges=[-2:0.1:5];   % edges to plot the lick rate histogram
    lick_trType=get_lickrate_byTrialType(trialData,trials,trialType,edges);
    plot_lickrate_byTrialType(lick_trType);
    
    %% plot response times
    valLabel='Response time (s)';    
    trialType={'go','left','right'};
    edges=[-1:0.05:2];
    respTime_trType=get_val_byTrialType(trialData.rt,trials,trialType,edges,valLabel);
    plot_val_byTrialType(respTime_trType);
    
    print(gcf,'-dpng','rt_byTrialType');    %png format
    saveas(gcf, 'rt_byTrialType', 'fig');

    %% plot ITI
    valLabel='Inter-trial interval (s)';    
    trialType={'go','reward','noreward'};
    edges=[0:0.25:20];
    iti_dur = [trialData.cueTimes(2:end)-trialData.outcomeTimes(1:end-1); NaN];  %ill-defined ITI for last trial
    iti_trType=get_val_byTrialType(iti_dur,trials,trialType,edges,valLabel);
    plot_val_byTrialType(iti_trType);
    
    print(gcf,'-dpng','iti_byTrialType');    %png format
    saveas(gcf, 'iti_byTrialType', 'fig');

    disp(['Number of ITI > 10 s = ' int2str(sum(iti_dur>10))]);
    disp(['Number of ITI > 20 s = ' int2str(sum(iti_dur>20))]);
    disp(['Number of ITI > 30 s = ' int2str(sum(iti_dur>30))]);
    
    %% fit with logistic regression
    num_regressor=15;
    x=1;    %analyze player 1
    
    %regressors = rewarded choice, unrewarded choice
    [lregRCUC_output, ~, ~, ~]=logreg_RCUC(stats,1,num_regressor);
    plot_logreg(lregRCUC_output,tlabel);
    print(gcf,'-dpng','logregRCUC');    %png format
    saveas(gcf, 'logregRCUC', 'fig');
    
    %regressors = choice, choice x reward
    [lregCRInt_output, ~, ~, ~]=logreg_CRInt(stats,1,num_regressor);
    plot_logreg(lregCRInt_output,tlabel);
    print(gcf,'-dpng','logregCRInt');    %png format
    saveas(gcf, 'logregCRInt', 'fig');   
    
    %% fit to WSLS and Q-learning algorithms
    
    fun = 'WSLSfun';
    initpar=0.5; % initial [prob_WSLS]
    lb=0;
    ub=1;
    [fitpar.WSLS, ~, bic.WSLS, nlike.WSLS]=fit_fun(stats,fun,initpar,1,lb,ub);
    
    fun = 'Q_RPEfun';
    initpar=[0.5 10]; % initial [alpha beta]
    [fitpar.Q_RPE, ~, bic.Q_RPE, nlike.Q_RPE]=fit_fun(stats,fun,initpar,1);
    
    fun = 'DFQfun';
    initpar=[0.1 0 0.8 0]; % initial [alpha1 alpha2 kappa1 kappa2]
    [fitpar.DFQ, ~, bic.DFQ, nlike.DFQ]=fit_fun(stats,fun,initpar,1);

    fun = 'FQfun';
    initpar=[0.1 0.8 0]; % initial [alpha1 kappa1 kappa2]
    [fitpar.FQ, ~, bic.FQ, nlike.FQ]=fit_fun(stats,fun,initpar,1);

    fun = 'Qfun';
    initpar=[0.1 0.8 0]; % initial [alpha1 kappa1 kappa2]
    [fitpar.Q, ~, bic.Q, nlike.Q]=fit_fun(stats,fun,initpar,1);
    
    [~, ~, bic.logregCRInt2, nlike.logregCRInt2]=logreg_CRInt(stats,1,2);
    [~, ~, bic.logregCRInt5, nlike.logregCRInt5]=logreg_CRInt(stats,1,5);
    [~, ~, bic.logregCRInt10, nlike.logregCRInt10]=logreg_CRInt(stats,1,10);

    fitpar
    bic
    nlike
    
    %% compare simulated choice behavior to actual choice behavior
    % without cross-validation, this plot may suffer from over-fitting

    player1.label='algo_DFQ';
    player1.params.a1=fitpar.FQ(1);    % learning rate (also = 1 minus the forgetting rate)
    player1.params.a2=fitpar.FQ(1);    % related to forgetting rate for action not chosen
    player1.params.k1=fitpar.FQ(2);    % strength of reinforcement by reward
    player1.params.k2=fitpar.FQ(3);    % strength of aversion from no-reward outcome
    
    stats_sim=predictAgent(player1,stats);
    
    plot_session_PrL(stats,stats_sim,sessionData.nTrials,tlabel);

    %%
    save(fullfile(savematpath,'beh.mat'),...
        'trialData','sessionData','trials',...
        'sw_output','sw_hrside_output','bl_output',...
        'lick_trType','iti_trType','respTime_trType',...
        'lregRCUC_output','lregCRInt_output',...
        'fitpar','bic','nlike','stats');

    close all;
    clearvars -except i dirs expData;
end

% plays sound when done
load train;
sound(y,Fs);
