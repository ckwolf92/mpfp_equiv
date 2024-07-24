%% TWO-TYPE FISCAL POLICY TRANSITION PATHS
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
model = '/twotype_lab';
task = '/ge';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path model '/_aux']);
addpath([path model '/_inputs/_results']);
addpath([path model task '/_results']);
cd([path model task]);

%% PREPARATIONS

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

global varphi delta chi_B chi_H epsilon_p phi_p kappa_p ...
     labtax rho_tr phi_pi phi_y phi_dy rho_dr

load param_agg

global beta gamma

load param_households

global C_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS B_SS D_SS Z_SS r_b_SS

load SS

global eta alpha mu

load aux

global C_w B_w L_w C_tau B_tau L_tau C_ib B_ib L_ib C_pi B_pi L_pi C_d B_d L_d C_b B_b L_b

load jac_matrices

%----------------------------------------------------------------
% Time Horizon
%----------------------------------------------------------------

global T

T = size(C_w,1);

%----------------------------------------------------------------
% NPV-maker
%----------------------------------------------------------------

r_NPV = zeros(T,1);
for t = 1:T
    r_NPV(t) = (1/(1+r_b_SS))^(t-1);
end

%% POLICY MATRICES

step = 1;

%----------------------------------------------------------------
% Monetary Policy Rule
%----------------------------------------------------------------

global dib_dpi dib_dy dib_dm

dib_dpi = NaN(T,T);

for i_deriv = 1:T
    pi_seq_deriv   = zeros(T,1);
    y_seq_deriv    = zeros(T,1);
    m_seq_deriv    = zeros(T,1);
    pi_seq_deriv(i_deriv,1) = pi_seq_deriv(i_deriv,1) + step;
    ib_up = mp_fn(pi_seq_deriv,y_seq_deriv,m_seq_deriv);
    dib_dpi(:,i_deriv) = ib_up/step;
end

dib_dy = NaN(T,T);

for i_deriv = 1:T
    pi_seq_deriv   = zeros(T,1);
    y_seq_deriv    = zeros(T,1);
    m_seq_deriv    = zeros(T,1);
    y_seq_deriv(i_deriv,1) = y_seq_deriv(i_deriv,1) + step;
    ib_up = mp_fn(pi_seq_deriv,y_seq_deriv,m_seq_deriv);
    dib_dy(:,i_deriv) = ib_up/step;
end

dib_dm = NaN(T,T);

for i_deriv = 1:T
    pi_seq_deriv   = zeros(T,1);
    y_seq_deriv    = zeros(T,1);
    m_seq_deriv    = zeros(T,1);
    m_seq_deriv(i_deriv,1) = m_seq_deriv(i_deriv,1) + step;
    ib_up = mp_fn(pi_seq_deriv,y_seq_deriv,m_seq_deriv);
    dib_dm(:,i_deriv) = ib_up/step;
end

%----------------------------------------------------------------
% Government Financing Rule
%----------------------------------------------------------------

global dtaue_dtaux dbg_dtaux dtaue_dw dbg_dw dtaue_dl dbg_dl dtaue_dpi dbg_dpi dtaue_dib dbg_dib

dtaue_dtaux = NaN(T,T);
dbg_dtaux   = NaN(T,T);

for i_deriv = 1:T
    taux_seq_deriv = zeros(T,1);
    ib_seq_deriv   = zeros(T,1);
    pi_seq_deriv   = zeros(T,1);
    w_seq_deriv    = zeros(T,1);
    l_seq_deriv    = zeros(T,1);
    taux_seq_deriv(i_deriv,1) = taux_seq_deriv(i_deriv,1) + step;
    [taue_up,bg_up] = gov_bc_fn(taux_seq_deriv,ib_seq_deriv,pi_seq_deriv,w_seq_deriv,l_seq_deriv);
    dtaue_dtaux(:,i_deriv) = taue_up/step;
    dbg_dtaux(:,i_deriv) = bg_up/step;
end

dtaue_dw  = NaN(T,T);
dbg_dw    = NaN(T,T);

for i_deriv = 1:T
    taux_seq_deriv = zeros(T,1);
    ib_seq_deriv   = zeros(T,1);
    pi_seq_deriv   = zeros(T,1);
    w_seq_deriv    = zeros(T,1);
    l_seq_deriv    = zeros(T,1);
    w_seq_deriv(i_deriv,1) = w_seq_deriv(i_deriv,1) + step;
    [taue_up,bg_up] = gov_bc_fn(taux_seq_deriv,ib_seq_deriv,pi_seq_deriv,w_seq_deriv,l_seq_deriv);
    dtaue_dw(:,i_deriv) = taue_up/step;
    dbg_dw(:,i_deriv) = bg_up/step;
end

dtaue_dl = dtaue_dw;
dbg_dl   = dbg_dw;

dtaue_dpi = NaN(T,T);
dbg_dpi   = NaN(T,T);

for i_deriv = 1:T
    taux_seq_deriv = zeros(T,1);
    ib_seq_deriv   = zeros(T,1);
    pi_seq_deriv   = zeros(T,1);
    w_seq_deriv    = zeros(T,1);
    l_seq_deriv    = zeros(T,1);
    pi_seq_deriv(i_deriv,1) = pi_seq_deriv(i_deriv,1) + step;
    [taue_up,bg_up] = gov_bc_fn(taux_seq_deriv,ib_seq_deriv,pi_seq_deriv,w_seq_deriv,l_seq_deriv);
    dtaue_dpi(:,i_deriv) = taue_up/step;
    dbg_dpi(:,i_deriv) = bg_up/step;
end

dtaue_dib = NaN(T,T);
dbg_dib   = NaN(T,T);

for i_deriv = 1:T
    taux_seq_deriv = zeros(T,1);
    ib_seq_deriv   = zeros(T,1);
    pi_seq_deriv   = zeros(T,1);
    w_seq_deriv    = zeros(T,1);
    l_seq_deriv    = zeros(T,1);
    ib_seq_deriv(i_deriv,1) = ib_seq_deriv(i_deriv,1) + step;
    [taue_up,bg_up] = gov_bc_fn(taux_seq_deriv,ib_seq_deriv,pi_seq_deriv,w_seq_deriv,l_seq_deriv);
    dtaue_dib(:,i_deriv) = taue_up/step;
    dbg_dib(:,i_deriv) = bg_up/step;
end

%% MODEL SOLUTION FOR TRANSFER POLICY

%----------------------------------------------------------------
% Shocks
%----------------------------------------------------------------

global zeta_seq m_seq taux_seq

zeta_seq = zeros(T,1);
m_seq    = zeros(T,1);

load results_mp

demand_target = (C_ib + C_tau * dtaue_dib) * dib_dm * m_seq_mp;

taux_seq = (C_tau + C_tau * dtaue_dtaux)^(-1) * demand_target;

demand_check = (C_tau + C_tau * dtaue_dtaux) * taux_seq;

%----------------------------------------------------------------
% Equilibrium Solution
%----------------------------------------------------------------

% initial guess

guess_seq = zeros(2*T,1);
excess_demand_init = excess_demand_fn(guess_seq);

% updating matrix

step = 10^(-3);
A_upd = NaN(2*T,2*T);

for i_deriv = 1:2*T
    guess_seq_deriv = zeros(2*T,1);
    guess_seq_deriv(i_deriv,1) = guess_seq_deriv(i_deriv,1) + step;
    excess_demand_up = excess_demand_fn(guess_seq_deriv);
    A_upd(:,i_deriv) = (excess_demand_up - excess_demand_init)/step;
end

% solution

x_sol = -A_upd\excess_demand_init;

pi_seq  = x_sol(1:T,1);
l_seq   = x_sol(T+1:2*T,1);

get_aggregates

%----------------------------------------------------------------
% Collect Results
%----------------------------------------------------------------

c_seq_fp      = c_seq;
y_seq_fp      = y_seq;
pi_seq_fp     = pi_seq;
b_seq_fp      = b_seq;
ib_seq_fp     = ib_seq;

%% PLOT RESULTS

cd([path model task '/_results']);

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

% IRF horizon

max_hor = 21;

% consumption impulse response

figure(1)
pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:max_hor-1,c_seq_mp(1:max_hor),'linewidth',5,'linestyle','-','color',settings.colors.grey)
hold on
plot(0:1:max_hor-1,c_seq_fp(1:max_hor),'linewidth',5,'linestyle',':','color',settings.colors.navyblue)
hold on
set(gcf,'color','w')
xlabel('Horizon','interpreter','latex','FontSize',18)
ylabel('\%','interpreter','latex','FontSize',18)
legend({'Interest Rate Cut','Stimulus Checks'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto');
print('figure_c5_2','-dpng');

cd([path model task]);