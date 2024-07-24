%% HANK FISCAL POLICY TRANSITION PATHS
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
model = '/hank/baseline';
task = '/ineq';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path model '/_aux']);
addpath([path model '/_inputs/_results']);
addpath([path model task '/_results']);
cd([path model task]);

%% PREPARATIONS

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

global beta beta_hat gamma probdeath varphi epsilon_w phi_w kappa_w wealth_0_pos ...
     epsilon_p phi_p kappa_p ...
     labtax TransY_ratio BY_ratio rho_tr phi_pi phi_y phi_dy rho_dr

load param_agg

global a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap borr_wedge

load param_households

global C_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS r_b_SS B_SS D_SS Z_SS ...
    C_sd_SS wealth_pctls C_pctl_SS liqwealth_indic_SS lambda_pctl_SS n_pctls ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

load SS

global grid_a spliorder states states_div states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

load aux

global C_w B_w C_l B_l C_tau B_tau C_ib B_ib C_pi B_pi C_d B_d C_b B_b

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

%% CONSTRUCT EQUIVALENT TRANSFER POLICY

%----------------------------------------------------------------
% Original Monetary Policy
%----------------------------------------------------------------

load results_mp

m_seq = m_seq_mp;

ib_seq_mp_PE   = dib_dm * m_seq;
taue_seq_mp_PE = dtaue_dib * ib_seq_mp_PE;

c_seq_mp_PE = C_ib * ib_seq_mp_PE + C_tau * taue_seq_mp_PE;

% implied inequality path

ib_seq    = ib_seq_mp_PE;
pi_seq    = zeros(T,1);
w_seq     = zeros(T,1);
l_seq     = zeros(T,1);
trans_seq = taue_seq_mp_PE;
d_seq     = zeros(T,1);
zeta_seq  = zeros(T,1);

solve_hh_problem

scale_approx = c_seq_mp_PE(1)/c_seq(1); % slight re-scaling due to finite differences

c_seq_mp_PE_approx    = scale_approx * c_seq;
c_sd_seq_mp_PE_approx = scale_approx * c_sd_seq;
c_pctl_seq_mp_PE      = scale_approx * c_pctl_seq;

%----------------------------------------------------------------
% Recover Equivalent Path
%----------------------------------------------------------------

tau_equiv_seq = C_tau^(-1) * c_seq_mp_PE;

% ensure 0 NPV (slight approximation error due to finite differences)

discount = NaN(T,1);
for t = 1:T
    discount(t) = 1/(1 + r_b_SS)^(t-1);
end
offset = - 1/sum(discount) * sum(discount .* tau_equiv_seq);
tau_equiv_seq = tau_equiv_seq - 1/sum(discount) * sum(discount .* tau_equiv_seq);

% translate to $ terms

tau_equiv_seq_dollar = tau_equiv_seq * Trans_SS/C_SS * 15000 / 100;

% implied inequality path

ib_seq    = zeros(T,1);
pi_seq    = zeros(T,1);
w_seq     = zeros(T,1);
l_seq     = zeros(T,1);
trans_seq = tau_equiv_seq;
d_seq     = zeros(T,1);
zeta_seq  = zeros(T,1);

solve_hh_problem

scale_approx = c_seq_mp_PE(1)/c_seq(1); % slight re-scaling due to finite differences

c_seq_fp_PE_approx    = scale_approx * c_seq;
c_sd_seq_fp_PE_approx = scale_approx * c_sd_seq;
c_pctl_seq_fp_PE      = scale_approx * c_pctl_seq;

%----------------------------------------------------------------
% Construct Inequality for Equivalent Policies
%----------------------------------------------------------------

c_pctl_seq_fp = c_pctl_seq_fp_PE + (c_pctl_seq_mp - c_pctl_seq_mp_PE);

%% PLOT RESULTS

%----------------------------------------------------------------
% Re-Scale Results
%----------------------------------------------------------------

scale = 1/c_seq_mp_PE(1);

c_seq_mp      = scale * c_seq_mp;
pi_seq_mp     = scale * pi_seq_mp;
ib_seq_mp_PE  = scale * ib_seq_mp_PE;
c_seq_mp_PE   = scale * c_seq_mp_PE;

c_seq_mp_PE_dollar          = c_seq_mp_PE * 15000/100;
tau_equiv_seq_dollar        = scale * tau_equiv_seq_dollar;

c_pctl_seq_mp_PE = scale * c_pctl_seq_mp_PE;
c_pctl_seq_mp    = scale * c_pctl_seq_mp;

c_pctl_seq_fp_PE = scale * c_pctl_seq_fp_PE;
c_pctl_seq_fp    = scale * c_pctl_seq_fp;

%----------------------------------------------------------------
% Color Preparation
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.lpurple = 0.25 * [160/255 32/255 240/255] + 0.75 * [1 1 1];
settings.colors.lorange = 0.25 * [204/255 102/255 0/255] + 0.75 * [1 1 1];

n = 200;

clear cmap
cmap(1,:) = [204/255,102/255,0/255];
cmap(2,:) = [1 1 1];
cmap(3,:) = [160/255,160/255,160/255];
[X,Y] = meshgrid([1:3],[1:50]);
cmap = interp2(X([1,25,50],:),Y([1,25,50],:),cmap,X,Y);

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

IRF_plot = 20;

%----------------------------------------------------------------
% Plot
%----------------------------------------------------------------

cd([path model task '/_results']);

plotwidth = 0.4;
gapsize = 0.075;
gapsize_edges = (1-2*plotwidth-1*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

figure(1)

subplot(1,2,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(wealth_pctls,[c_pctl_seq_mp_PE(1,:),c_pctl_seq_mp_PE(1,end)],'linewidth',5,'linestyle','-','color',settings.colors.green)
hold on
jbfill(wealth_pctls,[c_pctl_seq_mp_PE(1,:),c_pctl_seq_mp_PE(1,end)],...
    [c_pctl_seq_mp(1,:),c_pctl_seq_mp(1,end)],settings.colors.lpurple,settings.colors.lpurple,0,1);
hold on
plot(wealth_pctls,[c_pctl_seq_mp_PE(1,:),c_pctl_seq_mp_PE(1,end)],'linewidth',5,'linestyle','-','color',settings.colors.green)
hold on
set(gcf,'color','w')
title('Interest Rate Policy','interpreter','latex','fontsize',26)
xlabel('Wealth Percentile','interpreter','latex','FontSize',22)
ylabel('\% deviation','interpreter','latex','FontSize',22)
legend({'Direct','+ GE'},'Location','Northeast','fontsize',22,'interpreter','latex')
xlim([0 1])
ylim([0 2])
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
plot(wealth_pctls,[c_pctl_seq_fp_PE(1,:),c_pctl_seq_fp_PE(1,end)],'linewidth',5,'linestyle','-','color',settings.colors.navyblue)
hold on
jbfill(wealth_pctls,[c_pctl_seq_fp_PE(1,:),c_pctl_seq_fp_PE(1,end)],...
    [c_pctl_seq_fp(1,:),c_pctl_seq_fp(1,end)],settings.colors.lpurple,settings.colors.lpurple,0,1);
hold on
plot(wealth_pctls,[c_pctl_seq_fp_PE(1,:),c_pctl_seq_fp_PE(1,end)],'linewidth',5,'linestyle','-','color',settings.colors.navyblue)
hold on
set(gcf,'color','w')
title('Transfer Policy','interpreter','latex','fontsize',26)
xlabel('Wealth Percentile','interpreter','latex','FontSize',22)
legend({'Direct','+ GE'},'Location','Northeast','fontsize',22,'interpreter','latex')
xlim([0 1])
ylim([0 12])
grid on
hold off
pos = get(gcf, 'Position');

set(gcf, 'Position', [pos(1) pos(2) 1.9*pos(3) 1.15*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_c6','-depsc');

cd([path model task]);