%% Simulate a matching pennies game

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

n=5000;       % number of trials to simulate
n_plot=500;    % number of trials to plot

m = 13;

switch m
    case 1  
        player1.label='algo_stochastic';
        player1.params.p=0.6;          % probability of choosing left (e.g., 0.6, because there is some tolerance for deviation from 0.5)
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['Stochastic Prob(L)=' num2str(player1.params.p) ' versus Algorithm 2, n=' num2str(n) ' trials'];

    case 2  
        player1.label='algo_WSLS';
        player1.params.p=1;          
        player2.label='algo1';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['WSLS versus Algorithm 1, n=' num2str(n) ' trials'];
        
    case 3
        player1.label='algo_WSLS';
        player1.params.p=1;          
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['WSLS versus Algorithm 2, n=' num2str(n) ' trials'];
        
    case 4  %logistic regression; regressors=choice and choice x reward
        player1.label='algo_logreg_CRInt';
        player1.params.bias=0;          
        player1.params.Ch=[-0.02 0.16 0.13 0.14 0.12 0.11 0.09 0.08 0.09 0];         
        player1.params.RC=[0.22 0.08 0 0 0 0 0 0 0 0];                  
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['Logistic regression versus Algorithm 2, n=' num2str(n) ' trials'];

    case 5  %logistic regression; only the choice regressor
        %logit equation suggests these coefficient alone can lead to maximal P(Left) = 72% chance if all relevant prior choices were Left
        %steady state is 64% for one side, logit(0.64)=coeff*0.64
        player1.label='algo_logreg_CRInt';
        player1.params.bias=0;          
        player1.params.Ch=[0 0.16 0.13 0.14 0.12 0.11 0.09 0.08 0.09 0];         
        player1.params.RC=[0 0 0 0 0 0 0 0 0 0];                  
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['Logistic regression versus Algorithm 2, n=' num2str(n) ' trials'];

    case 6
        player1.label='algo_Q_RPE';
        player1.params.a=0.1;    % learning rate (e.g., 0.1 or lower to even 0)
        player1.params.b=5;      % inverse temperature (e.g., 1)
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['Q-learning w/ RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ...
            ' versus Algorithm 2, n=' num2str(n) ' trials'];

    case 7
        player1.label='algo_DFQ';
            %differential, forgetting, and Q-learning
            %forgetting and Q-learning (if a2==a1)
            %Q-learning (if a2==0)
        player1.params.a1=0.06;    % learning rate (also = 1 minus the forgetting rate)
        player1.params.a2=0.06;    % related to forgetting rate for action not chosen
        player1.params.k1=2.01;    % strength of reinforcement by reward
        player1.params.k2=-1.04;   % strength of aversion from no-reward outcome

        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['DF-Q a1=' num2str(player1.params.a1) ',a2=' num2str(player1.params.a2) ...
            ',k1=' num2str(player1.params.k1) ',k2=' num2str(player1.params.k2) ...
            ' versus Algorithm 2, n=' num2str(n) ' trials'];
    
    case 8
        player1.label='algo_DFQ_withbeta';
        %differential, forgetting, and Q-learning
        %forgetting and Q-learning (if a2==a1)
        %Q-learning (if a2==0)
        player1.params.a1=0.06;    % learning rate (also = 1 minus the forgetting rate)
        player1.params.a2=0.06;    % related to forgetting rate for action not chosen
        player1.params.k1=2.01;    % strength of reinforcement by reward
        player1.params.k2=-1.04;   % strength of aversion from no-reward outcome
        player1.params.beta = 2;
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['DF-Q a1=' num2str(player1.params.a1) ',a2=' num2str(player1.params.a2) ...
            ',k1=' num2str(player1.params.k1) ',k2=' num2str(player1.params.k2) ...
            ' versus Algorithm 2, n=' num2str(n) ' trials'];
    
    case 9
        player1.label='algo_DFQ_withbetabias';
        %differential, forgetting, and Q-learning
        %forgetting and Q-learning (if a2==a1)
        %Q-learning (if a2==0)
        player1.params.a1=0.06;    % learning rate (also = 1 minus the forgetting rate)
        player1.params.a2=0.06;    % related to forgetting rate for action not chosen
        player1.params.k1=2.01;    % strength of reinforcement by reward
        player1.params.k2=-1.04;   % strength of aversion from no-reward outcome
        player1.params.beta = 2;
        player1.params.deltaQ = 0.1;
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['DF-Q a1=' num2str(player1.params.a1) ',a2=' num2str(player1.params.a2) ...
            ',k1=' num2str(player1.params.k1) ',k2=' num2str(player1.params.k2) ...
            ' versus Algorithm 2, n=' num2str(n) ' trials'];
    case 10
        player1.label='algo_DFQ_withbetaCA';
        %differential, forgetting, and Q-learning
        %forgetting and Q-learning (if a2==a1)
        %Q-learning (if a2==0)
        player1.params.a1=0.06;    % learning rate (also = 1 minus the forgetting rate)
        player1.params.a2=0.06;    % related to forgetting rate for action not chosen
        player1.params.k1=2.01;    % strength of reinforcement by reward
        player1.params.k2=-1.04;   % strength of aversion from no-reward outcome
        player1.params.beta = 0.5;
        player1.params.tau = 0.5;
        player1.params.phi = 0.6;
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['DF-Q a1=' num2str(player1.params.a1) ',a2=' num2str(player1.params.a2) ...
            ',k1=' num2str(player1.params.k1) ',k2=' num2str(player1.params.k2) ...
            ' versus Algorithm 2, n=' num2str(n) ' trials'];
    case 11
        player1.label='algo_DFQ';
        player1.params.a1=0.06;    % learning rate (also = 1 minus the forgetting rate)
        player1.params.a2=0.06;    % related to forgetting rate for action not chosen
        player1.params.k1=2.01;    % strength of reinforcement by reward
        player1.params.k2=-1.04;   % strength of aversion from no-reward outcome
        
        player1.frac_opto=0.5;    % fraction of trial with optogenetic stim.
        player1.opto.k1=0;        % do not update for positive reinforcement
                
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['DF-Q (opto) a1=' num2str(player1.params.a1) ',a2=' num2str(player1.params.a2) ...
            ',k1=' num2str(player1.params.k1) ',k2=' num2str(player1.params.k2) ...
            ',fracopto=' num2str(player1.frac_opto) ',a1opto=' num2str(player1.opto.a1) ...           
            ' versus Algorithm 2, n=' num2str(n) ' trials'];
    case 12 
        player1.label='algo_DFQ';
        player1.params.a1=0.06;    % learning rate (also = 1 minus the forgetting rate)
        player1.params.a2=0.06;    % related to forgetting rate for action not chosen
        player1.params.k1=2.01;    % strength of reinforcement by reward
        player1.params.k2=-1.04;   % strength of aversion from no-reward outcom
                
        player2.label='algo_DFQ';
        player2.params.a1=0.1;   
        player2.params.a2=0.1;   
        player2.params.k1=2.01;    % strength of reinforcement by reward
        player2.params.k2=-1.04;   % strength of aversion from no-reward outcome
        
        tlabel=['DF-Q a1=' num2str(player1.params.a1) ',a2=' num2str(player1.params.a2) ...
            ',k1=' num2str(player1.params.k1) ',k2=' num2str(player1.params.k2) ...
            ' versus RL-DFQ, n=' num2str(n) ' trials'];
        
    case 13
        % set the payoff matrix to be
        
        % (2, -2)     (0, 0)
        % (1, -1)     (2, -2)  (mouse, computer)
        % nash equilibrium for the mouse is P_L = 2/3
        % nash equilibrium for the computer is P_L 1/3
        player1.label='algo_FQ_RPE';
        player1.params.a=0.5;    % learning rate (also = 1 minus the forgetting rate)
        player1.params.b = 2;
       
        %player1.bias = 1;  % if it is a biased payoff matrix
        
        player2.label='algo2_biased';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['F-Q-RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ... 
            ' versus algorithm2, n=' num2str(n) ' trials'];
end

%% simulate the game

disp('--- Simulating MP game ---');
disp(tlabel);
disp('--------------------------');
stats=simPennies(player1,player2,n);

%% plot simulated choice behavior - whole sessions
cd(savesimfigpath);

plot_session_game(stats,n_plot,tlabel);

%% plot latent learning parameters for player 1

x = 1;  %player 1
if strcmp(player1.label,'algo_Q_RPE') || strcmp(player1.label,'algo_DFQ')
    plot_session_qparam(stats,x,n_plot);
end

%% plot logistic regression

x = 1;  %player 1
num_regressor=15;

if isfield(player1,'frac_opto')
    [lreg_output, ~, ~, ~]=logreg_RCCh_opto(stats,x,num_regressor);
    plot_logreg_opto(lreg_output,tlabel);
    print(gcf,'-dpng','logregRCCh_opto');    %png format
    saveas(gcf, 'logregRCCh_opto', 'fig');   
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
[fitpar.Q_RPE, ~, bic.Q_RPE, nlike.Q_RPE]=fit_fun(stats,fun,initpar,1);

fun = 'DFQfun';
initpar=[0.1 0 0.8 0]; % initial [alpha1 alpha2 kappa1 kappa2]
[fitpar.DFQ, ~, bic.DFQ, nlike.DFQ]=fit_fun(stats,fun,initpar,1);

fun = 'FQfun';
initpar=[0.1 0.8 0]; % initial [alpha1 kappa1 kappa2]
[fitpar.FQ, ~, bic.FQ, nlike.FQ]=fit_fun(stats,fun,initpar,1);

fun = 'FQfun_withbeta';
initpar=[0.1 0 1]; % initial [alpha1 kappa2 beta]
[fitpar.FQ_wB, ~, bic.FQ_wB, nlike.FQ_wB]=fit_fun(stats,fun,initpar,1);

fun = 'FQfun_withbetabias';
initpar=[0.1 0 1 0]; % initial [alpha1 kappa2 beta deltaQ]
[fitpar.FQ_wBB, ~, bic.FQ_wBB, nlike.FQ_wBB]=fit_fun(stats,fun,initpar,1);

fun = 'FQfun_withbetaCA';
initpar=[0.1 0 1 0 0]; % initial [alpha1 kappa2 beta tau phi]
[fitpar.FQ_wBC, ~, bic.FQ_wBC, nlike.FQ_wBC]=fit_fun(stats,fun,initpar,1);

fun = 'Qfun';
initpar=[0.1 0.8 0]; % initial [alpha1 kappa1 kappa2]
[fitpar.Q, ~, bic.Q, nlike.Q]=fit_fun(stats,fun,initpar,1);

[~, ~, bic.logregCRInt2, nlike.logregCRInt2]=logreg_CRInt(stats,1,2);
[~, ~, bic.logregCRInt5, nlike.logregCRInt5]=logreg_CRInt(stats,1,5);
[~, ~, bic.logregCRInt10, nlike.logregCRInt10]=logreg_CRInt(stats,1,10);

fitpar
bic

%% time how long the simulation took
toc