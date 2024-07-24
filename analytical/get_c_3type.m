%% GENERATE SPENDER-2-OLG HYBRID INPUTS FOR POLICY EQUIVALENCE ANALYSIS
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
experiment = '/analytical';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path experiment '/_aux']);
addpath([path experiment '/_inputs']);
cd([path experiment]);

%% PARAMETERS

%----------------------------------------------------------------
% Calibration Target
%----------------------------------------------------------------

global MPC_target

load fhn_results

MPC_target = cumC1(5:10)';

%----------------------------------------------------------------
% Preferences
%----------------------------------------------------------------

global beta omega omega_1 omega_2 sigma MPC_target mu_1 mu_2

beta  = 0.99^(0.25);

calib_fn = @(param) calib_fn_aux(param);

param_sol = fminsearch(calib_fn,[0.8 0.6 0.3 0.5]);

omega_1 = param_sol(1);
omega_2 = param_sol(2);
mu_1    = param_sol(3);
mu_2    = param_sol(4);

sigma   = 1;

%----------------------------------------------------------------
% Steady-State Objects
%----------------------------------------------------------------

global tau_y Y_SS C_SS r_SS D_SS D_OLG_SS

tau_y    = 1/3;
Y_SS     = 1;
C_SS     = Y_SS;
r_SS     = 1/beta - 1;
D_SS     = 1.04;
D_OLG_SS = 1/(mu_1 + mu_2) * D_SS;

%% SETTINGS

global T

T    = 500;
step = 1;

exo.y_hat     = zeros(T,1);
exo.i_hat     = zeros(T,1);
exo.pi_hat    = zeros(T,1);
exo.zeta_hat  = zeros(T,1);

r_NPV = zeros(T,1);
for t = 1:T
    r_NPV(t) = (1/(1+r_SS))^(t-1);
end

%% DERIVATIVE MATRICES: TYPE 1

%----------------------------------------------------------------
% Auxiliary Matrix
%----------------------------------------------------------------

omega = omega_1;

A = zeros(2*T,2*T); % order (c,d) & (BC, EE)
for t = 1:T
    A(t,t) = 1;
    A(t,T+t) = 1;
    if t > 1
        A(t,T+t-1) = -1/beta;
    end
end
for t = T+1:2*T
    if t == T+1
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-T+1) = - beta * omega;
    elseif t < 2*T
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-T+1) = - beta * omega;
        A(t,t-1) = - (1-beta*omega) * (1-omega) * 1/beta;
    elseif t == 2*T
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-1) = 1;
    end
end

A_inv = A^(-1);

%----------------------------------------------------------------
% Baseline
%----------------------------------------------------------------

[c_base,d_base] = c_hybrid_fn(exo,A_inv);

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

C_y_1 = NaN(T,T);
D_y_1 = NaN(T,T);

for t = 1:T
    exo.y_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_y_1(:,t) = (c_shock-c_base)/step;
    D_y_1(:,t) = (d_shock-d_base)/step;
    exo.y_hat(t) = 0;
end

%----------------------------------------------------------------
% Transfers
%----------------------------------------------------------------

C_tau_1 = C_y_1;
D_tau_1 = D_y_1;

%----------------------------------------------------------------
% Nominal Interest Rates
%----------------------------------------------------------------

C_i_1 = NaN(T,T);
D_i_1 = NaN(T,T);

for t = 1:T
    exo.i_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_i_1(:,t) = (c_shock-c_base)/step;
    D_i_1(:,t) = (d_shock-d_base)/step;
    exo.i_hat(t) = 0;
end

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_pi_1 = NaN(T,T);
D_pi_1 = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_pi_1(:,t) = (c_shock-c_base)/step;
    D_pi_1(:,t) = (d_shock-d_base)/step;
    exo.pi_hat(t) = 0;
end

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_d_1 = NaN(T,T);
D_d_1 = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_d_1(:,t) = (c_shock-c_base)/step;
    D_d_1(:,t) = (d_shock-d_base)/step;
    exo.zeta_hat(t) = 0;
end

%% DERIVATIVE MATRICES: TYPE 2

%----------------------------------------------------------------
% Auxiliary Matrix
%----------------------------------------------------------------

omega = omega_2;

A = zeros(2*T,2*T); % order (c,d) & (BC, EE)
for t = 1:T
    A(t,t) = 1;
    A(t,T+t) = 1;
    if t > 1
        A(t,T+t-1) = -1/beta;
    end
end
for t = T+1:2*T
    if t == T+1
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-T+1) = - beta * omega;
    elseif t < 2*T
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-T+1) = - beta * omega;
        A(t,t-1) = - (1-beta*omega) * (1-omega) * 1/beta;
    elseif t == 2*T
        A(t,t-T) = 1 - omega * (1-beta*omega);
        A(t,t-1) = 1;
    end
end

A_inv = A^(-1);

%----------------------------------------------------------------
% Baseline
%----------------------------------------------------------------

[c_base,d_base] = c_hybrid_fn(exo,A_inv);

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

C_y_2 = NaN(T,T);
D_y_2 = NaN(T,T);

for t = 1:T
    exo.y_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_y_2(:,t) = (c_shock-c_base)/step;
    D_y_2(:,t) = (d_shock-d_base)/step;
    exo.y_hat(t) = 0;
end

%----------------------------------------------------------------
% Transfers
%----------------------------------------------------------------

C_tau_2 = C_y_2;
D_tau_2 = D_y_2;

%----------------------------------------------------------------
% Nominal Interest Rates
%----------------------------------------------------------------

C_i_2 = NaN(T,T);
D_i_2 = NaN(T,T);

for t = 1:T
    exo.i_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_i_2(:,t) = (c_shock-c_base)/step;
    D_i_2(:,t) = (d_shock-d_base)/step;
    exo.i_hat(t) = 0;
end

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_pi_2 = NaN(T,T);
D_pi_2 = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_pi_2(:,t) = (c_shock-c_base)/step;
    D_pi_2(:,t) = (d_shock-d_base)/step;
    exo.pi_hat(t) = 0;
end

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_d_2 = NaN(T,T);
D_d_2 = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_shock,d_shock] = c_hybrid_fn(exo,A_inv);
    C_d_2(:,t) = (c_shock-c_base)/step;
    D_d_2(:,t) = (d_shock-d_base)/step;
    exo.zeta_hat(t) = 0;
end

%% DERIVATIVE MATRICES: AGGREGATE

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

C_y = mu_1 * C_y_1 + mu_2 * C_y_2 + (1 - mu_1 - mu_2) * eye(T);
D_y = mu_1 * D_y_1 + mu_2 * D_y_2 + (1 - mu_1 - mu_2) * zeros(T,T);

%----------------------------------------------------------------
% Transfers
%----------------------------------------------------------------

C_tau = C_y;
D_tau = D_y;

%----------------------------------------------------------------
% Nominal Interest Rates
%----------------------------------------------------------------

C_i = mu_1 * C_i_1 + mu_2 * C_i_2 + (1 - mu_1 - mu_2) * zeros(T,T);
D_i = mu_1 * D_i_1 + mu_2 * D_i_2 + (1 - mu_1 - mu_2) * zeros(T,T);

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_pi = mu_1 * C_pi_1 + mu_2 * C_pi_2 + (1 - mu_1 - mu_2) * zeros(T,T);
D_pi = mu_1 * D_pi_1 + mu_2 * D_pi_2 + (1 - mu_1 - mu_2) * zeros(T,T);

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_d = mu_1 * C_d_1 + mu_2 * C_d_2 + (1 - mu_1 - mu_2) * zeros(T,T);
D_d = mu_1 * D_d_1 + mu_2 * D_d_2 + (1 - mu_1 - mu_2) * zeros(T,T);

%% COMPUTE CUMULATIVE SPENDING SHARES

cons_share = cumsum(r_NPV .* C_y(:,1));

cons_NPV = r_NPV .* C_y(:,1);
cons_share_annual = zeros(T/4,1);
for t = 1:T/4
    cons_share_annual(t) = sum(cons_NPV(1+(t-1)*4:t*4,1));
end
cons_share_annual = cumsum(cons_share_annual);

%% COMPUTE ACCURACY OF SUFFICIENT STATISTICS FORMULA

%----------------------------------------------------------------
% Get Sufficient Statistics
%----------------------------------------------------------------

C_tau_true = C_tau;
C_tau_true_inv = C_tau_true^(-1);

omega   = C_tau_true(1,1);
xitheta = C_tau_true(2,1)/omega;
theta   = (1 + r_SS) - xitheta * omega/(1-omega);
xi      = xitheta/theta;

psi = 1;

%----------------------------------------------------------------
% Get Approximation
%----------------------------------------------------------------

C_tau_pred = suffstats_fn(omega,theta,xi,psi,T);
C_tau_pred_inv = C_tau_pred^(-1);

%% PLOT RESULTS

cd([path experiment '/_results']);

% colors

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.lnavyblue = 0.25 * [0/255 0/255 50/255] + 0.75 * [1 1 1];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.lpurple = 0.25 * [160/255 32/255 240/255] + 0.75 * [1 1 1];

weights = [1 0.8 0.5 0.2];
settings.colors.orange_all = zeros(length(weights),3);
settings.colors.blue_all = zeros(length(weights),3);
for i = 1:4
    settings.colors.orange_all(i,:) = weights(i) * settings.colors.orange + (1-weights(i)) * [1 1 1];
    settings.colors.blue_all(i,:) = weights(i) * settings.colors.black + (1-weights(i)) * [1 1 1];
end

% show iMPC accuracy & suff stats accuracy

plotwidth = 0.385;
gapsize = 0.055;
gapsize_edges = (1-2*plotwidth-1*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

select_hor = [1 6 12 18];

max_hor = 21;

figure(1)

subplot(1,2,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
jbfill(-4:1:5,lowCI,uppCI,...
    0.9 * [1 1 1],0.9 * [1 1 1],0,1);
hold on
plot(-4:1:5,[zeros(4,1);cons_share_annual(1:6)],'linewidth',5,'linestyle','-','color',settings.colors.blue_all(2,:))
hold off
set(gcf,'color','w')
title('Cumulative MPCs','interpreter','latex','fontsize',21)
xlabel('Year','interpreter','latex','FontSize',18)
legend({'Data','Model'},'Location','Southeast','fontsize',18,'interpreter','latex')
ylim([-0.19 1])
xlim([-4 5])
grid on
hold off

subplot(1,2,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
for i = 1:length(select_hor)
    plot(0:1:max_hor-1,C_tau_true_inv(1:max_hor,select_hor(i)),'linewidth',5,'linestyle','-','color',settings.colors.blue_all(i,:))
    hold on
    plot(0:1:max_hor-1,C_tau_pred_inv(1:max_hor,select_hor(i)),'linewidth',5,'linestyle',':','color',settings.colors.orange_all(i,:))
    hold on
end
set(gcf,'color','w')
title('Entries of $\mathcal{C}_\tau^{-1}$','interpreter','latex','fontsize',21)
xlabel('Horizon','interpreter','latex','FontSize',18)
legend({'Exact Multi-Type','Sufficient Statistics'},'Location','Southeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) 1.2*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_c1','-depsc');

cd([path experiment]);