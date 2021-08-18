%% Simulate a biased matching pennies game

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

m = 1;

switch m
    case 1  
        player1.label='algo_stochastic';
        player1.params.p=0.7;          % probability of choosing left (e.g., 0.6, because there is some tolerance for deviation from 0.5)
        player2.label='algo2biased';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        player2.params.pLeft = 0.5;  %test deviation from this probability to choose left
        
        tlabel=['Stochastic Prob(L)=' num2str(player1.params.p) ' versus Algorithm 2, n=' num2str(n) ' trials'];

    case 2  %logistic regression; regressors=choice and choice x reward
        player1.label='algo_logreg_CRInt';
        player1.params.bias=0;          
        player1.params.Ch=[-0.02 0.16 0.13 0.14 0.12 0.11 0.09 0.08 0.09 0];         
        player1.params.RC=[0.22 0.08 0 0 0 0 0 0 0 0];                  
        player2.label='algo2biased';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        player2.params.pLeft = 0.75;  %test deviation from this probability to choose left
        
        tlabel=['Logistic regression versus Algorithm 2, n=' num2str(n) ' trials'];

    case 3
        player1.label='algo_FQ_RPE';
%             %differential, forgetting, and Q-learning
%             %forgetting and Q-learning (if a2==a1)
%             %Q-learning (if a2==0)
        player1.params.a=0.5;    % learning rate (also = 1 minus the forgetting rate)
% %        player1.params.k2=-1.04;   % strength of aversion from no-reward outcome
        player1.params.b=1;   % strength of aversion from no-reward outcome

        player2.label='algo2biased';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        player2.params.pLeft = 0.5;  %test deviation from this probability to choose left
        
        tlabel=['F-Q-RPE a=' num2str(player1.params.a) ',b=' num2str(player1.params.b) ... 
           ' versus algorithm2, n=' num2str(n) ' trials'];

end

%% simulate the game

disp('--- Simulating Biased MP game ---');
disp(tlabel);
disp('--------------------------');
payoff=[2 0;0 1];   %payoff matrix for player 1

% %parameter sweep 
% alist = 0.1:0.1:0.9;
% blist = 0:0.2:4;
% average_outcome = zeros(length(alist), length(blist));
% for ii = 1:length(alist)
%     player1.params.a=alist(ii);    % learning rate (also = 1 minus the forgetting rate)
%     parfor jj = 1:length(blist)
%         tic;
%         temp_player = struct();
%        temp_player.label='algo_FQ_RPE';
%             differential, forgetting, and Q-learning
%             forgetting and Q-learning (if a2==a1)
%             Q-learning (if a2==0)
%         tt = jj;
%         
%        player1.params.k2=-1.04;   % strength of aversion from no-reward outcome
%         temp_player.params.a=alist(ii); 
%         temp_player.params.b=blist(tt);   % strength of aversion from no-reward outcome
% 
%         
%         stats=simPenniesBiased(temp_player,player2,n,payoff);
%         average_outcome(ii,jj) = sum(stats.r(:,1))/n;
%         toc
%     end
% end
% 
% figure; surf(blist, alist, average_outcome)

stats=simPenniesBiased(player1,player2,n,payoff);
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