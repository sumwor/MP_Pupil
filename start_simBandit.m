%% Simulate a two-armed bandit task

clearvars;
close all;
tic;

%% setup path and plotting formats

value_setPathList;

setup_figprop;  %set up default figure plotting parameters

%%
savesimfigpath = fullfile(data_dir,'figs-sim');
if ~exist(savesimfigpath,'dir')
    mkdir(savesimfigpath);
end

%% set up opponents

n=100000;      % number of trials to simulate
n_plot=500;  % number of trials to plot

m = 6;

switch m
    
    case 1 %win-stay, lose-switch strategy, can get decent performance in bandit task
        player1.label='algo_WSLS';
        player1.params.p=0.9;    %probability of following WSLS         
        
        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];  
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};
        tlabel=['WSLS agent in two-armed bandit task, n=' num2str(n) ' trials'];

    case 2
        player1.label='algo_logreg_RCUC';
        player1.params.bias=0;          
        player1.params.RC=[2.5 1.1 0.6];          
        player1.params.UC=[-0.4 -0.1 0]; 
        
        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];  
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};
        tlabel=['Logistic regression agent'...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];

	case 3
        player1.label='algo_Q_RPE';
        player1.params.a=0.6;  % learning rate (e.g., 0.6)
        player1.params.b=5;    % inverse temperature (e.g., 4)
        
        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];  
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};
        tlabel=['Q-learning w/ RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];

    case 4
        player1.label='algo_DFQ';
            %differential, forgetting, and Q-learning
            %forgetting and Q-learning (if a2==a1)
            %Q-learning (if a2==0)
        player1.params.a1=0.5194;    % learning rate (also = 1 minus the forgetting rate)
        player1.params.a2=0.5194;    % related to forgetting rate for action not chosen
        player1.params.k1=3.6830;    % strength of reinforcement by reward
        player1.params.k2=-1.1264;   % strength of aversion from no-reward outcome

        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];  
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};

        tlabel=['DF-Q a1=' num2str(player1.params.a1) ',a2=' num2str(player1.params.a2) ...
            ',k1=' num2str(player1.params.k1) ',k2=' num2str(player1.params.k2) ...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];
    case 5
        player1.label='algo_associability';
            %Pearce-Hall associability

        player1.params.kappa=2;    % learning rate coefficient
        player1.params.eta=0.2;    % learning rate for updating associability
        player1.params.alpha0=0.2;    % initial value for associability

        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];  
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};

        tlabel=['Associability kappa=' num2str(player1.params.kappa) ',eta=' num2str(player1.params.eta) ...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];
    
    case 6
        player1.label='algo_Hybrid_FQBayesian';
        player1.params =[0 -1 1 0.3]; % a, k2, beta, w
%         player1.params.a=0;   
%         % player1.params.k1 = 2;
%         player1.params.k2 = -1;
%         player1.params.beta = 1;
%         player1.params.w = 0.3;

        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];  
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};

        tlabel=['hybrid kappa2=' num2str(player1.params(2)) ',w=' num2str(player1.params(4)) ...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];
    
    case 12
        player1.label='algo_logreg_RCUC';
        player1.params.bias=0;          
        player1.params.RC=[2.5 1.1 0.6];          
        player1.params.UC=[-0.4 -0.1 0]; 
        
        player1.frac_opto=0.5;    % fraction of trial with optogenetic stim.
        player1.opto.RC=[2.5 0 0];      % regression coefficients w/ optogenetic stim.

        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];  
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};
        tlabel=['Logistic regression agent, fracopto=' num2str(player1.frac_opto) ...           
            ' in two-armed bandit task, n=' num2str(n) ' trials'];
        
	case 13
        player1.label='algo_Q_RPE';
        player1.params.a=0.6;  % learning rate (e.g., 0.6)
        player1.params.b=5;    % inverse temperature (e.g., 4)
        
        player1.frac_opto=0.5;    % fraction of trial with optogenetic stim.
        player1.opto.a=1.2;       % regression coefficients w/ optogenetic stim.

        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];  
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};
        tlabel=['Q-learning w/ RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
            ',fracopto=' num2str(player1.frac_opto) ',aopto=' num2str(player1.opto.a) ...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];
    
    case 14  %Gershman, Cognition, 2017
        player1.label='algo_UCB';
        
        % par from the gershman paper (not suited)
        % player1.params.tau2=[10,10];  % learning rate (e.g., 0.6)
        % player1.params.sigma2=[100,100];    % inverse temperature (e.g., 4)
        
        %try different parameters
        player1.params.gamma = 0.5;  % for UCB uncertainty bonus (value?)
        player1.params.lambda = 0.25; % for UCB (value?);
        
        % assuming a beta prior
        % try different prios?
        player1.params.aL = 1;
        player1.params.bL = 1;
        player1.params.aR = 1;
        player1.params.bR = 1;
        %taskparams.p_pairs=[0.1 0.9; 0.1 0.5; 0.5 0.9; 0.5 0.5; 0.9 0.5; 0.5 0.1; 0.9 0.1]; 
        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};
        tlabel=['directed exploration: UCB ' ...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];
     
    case 15  %Gershman, Cognition, 2017
        player1.label='algo_Thompson';
        
        % par from the gershman paper (not suited)
        % player1.params.tau2=[10,10];  % learning rate (e.g., 0.6)
        % player1.params.sigma2=[100,100];    % inverse temperature (e.g., 4)
        
        % assuming a beta prior
        player1.params.aL = 1;
        player1.params.bL = 1;
        player1.params.aR = 1;
        player1.params.bR = 1;
        %taskparams.p_pairs=[0.1 0.9; 0.1 0.5; 0.5 0.9; 0.5 0.5; 0.9 0.5; 0.5 0.1; 0.9 0.1]; 
        taskparams.p_pairs=[0.7 0.1; 0.1 0.7];
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.7:0.1','0.1:0.7'};
        tlabel=['directed exploration: Thompson ' ...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];
    
    case 99  %Hamid et al., NN, 2016
        player1.label='algo_Q_RPE';
        player1.params.a=0.6;  % learning rate (e.g., 0.6)
        player1.params.b=4;    % inverse temperature (e.g., 4)
        
        taskparams.p_pairs=[0.1 0.9; 0.1 0.5; 0.5 0.9; 0.5 0.5; 0.9 0.5; 0.5 0.1; 0.9 0.1]; 
        taskparams.crit_hit=10;       % switching criterion: number of times picked the high-reward side
        taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
        rule_labels = {'0.1:0.9','0.1 0.5','0.5:0.9','0.5:0.5','0.9:0.5','0.5:0.1','0.9:0.1'};
        tlabel=['Q-learning w/ RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
            ' in two-armed bandit task, n=' num2str(n) ' trials'];

end

%% simulate the game

disp('--- Simulating bandit task ---');
disp(tlabel);
disp('------------------------------');
stats=simBandit(player1,taskparams,n);

%% plot simulated choice behavior - whole sessions
cd(savesimfigpath);

plot_session_task(stats,n_plot,tlabel);

%% plot simulated choice behavior - around switches

trials_back=10;  % set number of previous trials 

sw_output=choice_switch(stats,trials_back);
plot_switch(sw_output,tlabel,rule_labels);

sw_hrside_output=choice_switch_hrside(stats,trials_back);
plot_switch_hrside(sw_hrside_output,tlabel);

%% plot simulated choice behavior - during a block

trials_forw=10;  % set number of trials 

bl_output=choice_block(stats,trials_forw);
plot_block(bl_output,tlabel,rule_labels);

%% plot latent learning parameters for player 1

x = 1;  %player 1
if strcmp(player1.label,'algo_Q_RPE') || strcmp(player1.label,'algo_DFQ')
    plot_session_qparam(stats,x,n_plot); %t_lable
end

%% plot logistic regression

x = 1;  %player 1
num_regressor=15;
if isfield(player1,'frac_opto')
    [lreg_output, ~, ~, ~]=logreg_RCUC_opto(stats,x,num_regressor);
    plot_logreg_opto(lreg_output,tlabel);
    print(gcf,'-dpng','logregRCUC_opto');    %png format
    saveas(gcf, 'logregRCUC_opto', 'fig');   
else
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
end

%% fit to WSLS and Q-learning algorithms

fun = 'WSLSfun';
initpar=0.5; % initial [prob_WSLS]
lb=0;
ub=1;
[fitpar.WSLS, ~, bic.WSLS, nlike.WSLS]=fit_fun(stats,fun,initpar,1,lb,ub);

fun = 'Q_RPEfun';
initpar=[0.5 10]; % initial [alpha beta]
[fitpar.Q_RPE, ~, bic.Q_RPE, nlike.Q_RPE]=fit_fun(stats,fun,initpar,1)
% newstats.c= stats.c; newstats.r = stats.r;
% [fitpar.Q_RPE1, ~, bic.Q_RPE1, nlike.Q_RPE1]=fit_fun(newstats,fun,initpar,1);
fun = 'DFQfun';
initpar=[0.1 0 0.8 0]; % initial [alpha1 alpha2 kappa1 kappa2]
[fitpar.DFQ, ~, bic.DFQ, nlike.DFQ]=fit_fun(stats,fun,initpar,1);

fun = 'FQfun';
initpar=[0.1 0.8 0]; % initial [alpha1 kappa1 kappa2]
[fitpar.FQ, ~, bic.FQ, nlike.FQ]=fit_fun(stats,fun,initpar,1);

fun = 'Qfun';
initpar=[0.1 0.8 0]; % initial [alpha1 kappa1 kappa2]
[fitpar.Q, ~, bic.Q, nlike.Q]=fit_fun(stats,fun,initpar,1);

fun = 'FQfun_associability';
initpar=[0.7 0.1 0]; % initial [kappa eta alpha0]
[fitpar.associability, ~, bic.associability, nlike.associability]=fit_fun(stats,fun,initpar,1);

fun = 'funHybrid_FQBayesian';
initpar = [0 0 0 0];  % alpha, kappa2, beta w
[fitpar.hybrid, ~, bic.hybrid, nlike.hybrid]=fit_fun(stats,fun,initpar,1);

[~, ~, bic.logregCRInt2, nlike.logregCRInt2]=logreg_CRInt(stats,1,2);
[~, ~, bic.logregCRInt5, nlike.logregCRInt5]=logreg_CRInt(stats,1,5);
[~, ~, bic.logregCRInt10, nlike.logregCRInt10]=logreg_CRInt(stats,1,10);

fitpar
bic

%% plot the transition in different blocks
trans_block1 = [];
trans_block2 = [];
trans_block3 = [];
trans_block4 = [];

blockLength = [];
block = 1;
numSwitch = 0;
blockTrans = [];
rule = [];
for ii = 1:n-1
    if stats.rule(ii+1) == stats.rule(ii)
        block = block + 1;
    elseif stats.rule(ii+1) ~= stats.rule(ii)
        rule = [rule, stats.rule(ii)];
        blockTrans = [blockTrans, ii+1];
        blockLength = [blockLength, block];
        block = 0;
        numSwitch = numSwitch + 1;
    end
end

for ii = 1:length(blockTrans)
    if blockLength(ii)<= 25
        if rule(ii) == 1
            trans_block1 = [trans_block1; stats.c(blockTrans(ii)-10:blockTrans(ii)+10)'];
        else
            trans_block1 = [trans_block1; -stats.c(blockTrans(ii)-10:blockTrans(ii)+10)'];
        end
    elseif blockLength(ii) > 25 & blockLength(ii) <= 40
        if rule(ii) == 1
            trans_block2 = [trans_block2; stats.c(blockTrans(ii)-10:blockTrans(ii)+10)'];
        else
            trans_block2 = [trans_block2; -stats.c(blockTrans(ii)-10:blockTrans(ii)+10)'];
        end
    elseif blockLength(ii) > 40 & blockLength(ii) <= 60
        if rule(ii) == 1
            trans_block3 = [trans_block3; stats.c(blockTrans(ii)-10:blockTrans(ii)+10)'];
        else
            trans_block3 = [trans_block3; -stats.c(blockTrans(ii)-10:blockTrans(ii)+10)'];
        end
    elseif blockLength(ii) > 60 
        if rule(ii) == 1
            trans_block4 = [trans_block4; stats.c(blockTrans(ii)-10:blockTrans(ii)+10)'];
        else
            trans_block4 = [trans_block4; -stats.c(blockTrans(ii)-10:blockTrans(ii)+10)'];
        end
    end
        
end

% estimate choice probability
m1 = length(trans_block1); m2 = length(trans_block2); m3 = length(trans_block3); m4 = length(trans_block4);

p1 = sum(trans_block1 == -1)/m1;
p2 = sum(trans_block2 == -1)/m2;
p3 = sum(trans_block3 == -1)/m3;
p4 = sum(trans_block4 == -1)/m4;

figure; plot([-10:10],p1);
hold on; plot([-10:10],p2); hold on; plot([-10:10],p3); hold on; plot([-10:10],p4);
hold on; plot([0 0], [0 1],'--');
legend('10-25', '25-40','40-60','60-120');
%% fit the learning rate in a running window
window_step = 10;
window_length = 100;

learning_rate = [];
for ii = 1:window_step:n-window_length
    fun = 'Q_RPEfun';
    initpar=[0.5 10]; % initial [alpha beta]
    statsWin.c = stats.c(ii:ii+window_length);
    statsWin.r = stats.r(ii:ii+window_length);
    [fitpar.Q_RPE, ~, bic.Q_RPE, nlike.Q_RPE]=fit_fun(statsWin,fun,initpar,1);
    learning_rate = [learning_rate, fitpar.Q_RPE(1)];
end
figure;plot(learning_rate)

%% use python code to fit the baysian model
modPath = 'E:\labcode\optlearner-master\optlearner';
cd(modPath);
mod = py.importlib.import_module('optlearner');
py.importlib.import_module('seaborn');
learner = mod.VolatilityLearner();

%modify the stats.r
observation = zeros(1, length(stats.r));
observation(stats.c==-1) = stats.r(stats.c==-1);
observation(stats.c==1) = 1 - stats.r(stats.c==1);
learner.fit(observation);
learner.plot_history(stats.rewardprob(:,1)');
estimated_volatility = double(py.array.array('d', learner.I_hats));
estimated_probability = double(py.array.array('d', learner.p_hats));
estimated_k = double(py.array.array('d', learner.k_hats));

% get the joint distribution
shape_pI = double(py.array.array('d',learner.pI.shape));
joint_pI = double(py.array.array('d',py.numpy.nditer(learner.pI)));
joint_pI = reshape(joint_pI,shape_pI);
I_grid = cell(learner.I_dists);
%plot the choice, reward and estiamted volitility
trials = 1:500;
figure;
subplot(3,1,1)
plot(stats.rewardprob(trials,1),'r','LineWidth',2); hold on;
plot(estimated_probability(trials));
hold on; scatter(trials-trials(1), observation(trials));
subplot(3,1,2)
yyaxis left
plot(estimated_volatility(trials));
yyaxis right
hold on;plot([1:10:500],learning_rate(1:50))
subplot(3, 1, 3)
plot(estimated_k(trials));
% check the block volitility around block switch
numSwitch = 0;
estVol = [];
for ii = 1:n-1
    if stats.rule(ii+1) ~= stats.rule(ii)
        numSwitch = numSwitch + 1;
        estVol = [estVol; estimated_volatility(ii-10:ii+10)];
    end
end
figure
plot([-10:10], mean(estVol,1));
xlabel('Trials from the switch (0 is the first trial after switch)');
title('Volatility around block switch');

%% fit the regression for UCB and Thompson sampling agent
if m == 14 %(or 15 for Thompson)
    
    %get regression factors
    
    fun = 'BTHfun';
    initpar = [1, 1, 1];
    [fitpar.BTH, ~, bic.BTH, nlike.BTH, hess]=fit_Explore(stats,fun,initpar,1);
    err_UCB = sqrt(diag(inv(hess)));
    
    % in regression coefficient
    % w1 = 1/lambda; w2 = gamma/lambda
    % plot the result
    figure;
    bar(fitpar.BTH,'black');
    hold on;errorbar([1,2,3],fitpar.BTH,err_UCB,'.','LineWidth',2,'color','black');
    ylabel('Coefficient');
    xticklabels({'V','RU','V/TU'});
    title('Regression coefficient for UCB agent');
elseif m == 15
    
    %get regression factors
    
    fun = 'BTHfun';
    initpar = [1, 1, 1];
    [fitpar.BTH, ~, bic.BTH, nlike.BTH, hess]=fit_Explore_singleV(stats,fun,initpar,1);
    err_Thompson = sqrt(diag(inv(hess)));
    % plot the result
    figure;
    bar(fitpar.BTH', 'black');
    hold on;errorbar([1,2,3],fitpar.BTH,err_Thompson,'.','LineWidth',2,'color','black');
    ylabel('Coefficient');
    xticklabels({'V','RU','V/TU'});
    title('Regression coefficient for Thompson agent');
    
end

print(gcf,'-dpng','regression coefficient');    %png format
saveas(gcf, 'logregCRInt', 'fig');
%% time how long the simulation took
toc