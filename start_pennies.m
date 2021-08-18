%% Analyze performance of a matching pennies game

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

%% process data files
for i = 1:numel(expData)
    disp(['Processing file ' int2str(i) ' out of ' int2str(numel(expData)) '.']);
    disp(['Filename: ' expData(i).logfile ]);
    
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
    stats.playerlabel{2} = 'Algo2';
    
    %choice: left=-1; right=1; miss=NaN
    stats.c=nan(sessionData.nTrials,2);
    stats.c(trials.left,1)=-1;
    stats.c(trials.right,1)=1;
    stats.c(trials.compleft,2)=-1;
    stats.c(trials.compright,2)=1;
 
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
    
    plot_session_game(stats,sessionData.nTrials,tlabel);
    
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
    [lregRCUC_output, ~, ~, ~]=logreg_RCUC(stats,x,num_regressor);
    plot_logreg(lregRCUC_output,tlabel);
    print(gcf,'-dpng','logregRCUC');    %png format
    saveas(gcf, 'logregRCUC', 'fig');
    
    %regressors = choice, choice x reward (equivalent to computer's choice)
    [lregCRInt_output, ~, ~, ~]=logreg_CRInt(stats,x,num_regressor);
    plot_logreg(lregCRInt_output,tlabel);
    print(gcf,'-dpng','logregCRInt');    %png format
    saveas(gcf, 'logregCRInt', 'fig');
  
    %% fit to WSLS and Q-learning algorithms
    % normalized likelihood (Ito and Doya, PLoS Comp Biol, 2015)

    fun = 'WSLSfun';
    initpar=0.5; % initial [prob_WSLS]
    lb=0;
    ub=1;
    [fitpar.WSLS, ~, bic.WSLS, nlike.WSLS]=fit_fun(stats,fun,initpar,1,lb,ub);

    fun = 'Q_RPEfun';
    initpar=[0.5 10]; % initial [alpha beta]
    lb=[0 0];
    ub=[1 inf];
    [fitpar.Q_RPE, ~, bic.Q_RPE, nlike.Q_RPE]=fit_fun(stats,fun,initpar,1,lb,ub);

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
    
    %% compare simulated choice behavior to actual choice behavior
    % without cross-validation, this plot may suffer from over-fitting
    
    player1.label='algo_logreg_CRInt';
    player1.params.bias=lregCRInt_output.b_bias;
    player1.params.Ch=lregCRInt_output.b_coeff(:,1);
    player1.params.RC=lregCRInt_output.b_coeff(:,2);

    stats_sim=predictAgent(player1,stats);
    
    plot_session_PrL(stats,stats_sim,sessionData.nTrials,tlabel);
    
    %% calculate entropy
    
    % the possible choice combinations
    choiceBack = 3;
    combos = de2bi([0:2^choiceBack-1],choiceBack);
    combos = 2*(combos - 0.5);  %to make it -1 or 1
    
    % classify animal's choice sequence
    nTrial = size(stats.c,1);
    cumuoccur = zeros(2^choiceBack,nTrial);   %how many times each combo occurred, cumulatively
    for j = choiceBack:nTrial
        c = stats.c(j-choiceBack+1:j,1)';
        idx = ismember(combos,c,'rows');   %find if it matches one of the combos
        if sum(idx)==1
            cumuoccur(:,j) = cumuoccur(:,j-1) + idx;  % update the cumulative occurrence matrix
        else
            cumuoccur(:,j) = cumuoccur(:,j-1);
        end
    end
    
%     figure; hold on;
%     for k=1:2^choiceBack
%         plot(cumuoccur(k,:));
%     end
    
    p = cumuoccur(:,end)/sum(cumuoccur(:,end));
    entro = -1*sum(p.*log2(p));
    
    %%
    save(fullfile(savematpath,'beh.mat'),...
            'trialData','sessionData','trials',...
            'lregRCUC_output','lregCRInt_output',...
            'lick_trType','iti_trType','respTime_trType',...
            'fitpar','bic','nlike',...
            'entro',...
            'stats');

    close all;
    clearvars -except i dirs expData;
end

% plays sound when done
load train;
sound(y,Fs);
