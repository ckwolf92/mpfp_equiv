%% HANK MONETARY POLICY TRANSITION PATHS
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
model = '/hank/labor_variant';
task = '/ge';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path model '/_aux']);
addpath([path model '/_inputs/_results']);
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

global C_w Cs_w B_w C_l Cs_l B_l C_ib Cs_ib B_ib C_pi Cs_pi B_pi ...
    C_tau Cs_tau B_tau C_d Cs_d B_d C_b Cs_b B_b

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

%% MODEL SOLUTION FOR INTEREST RATE POLICY

%----------------------------------------------------------------
% Shocks
%----------------------------------------------------------------

global zeta_seq m_seq taux_seq

zeta_seq = zeros(T,1);
m_seq    = zeros(T,1);
taux_seq = zeros(T,1);

m_seq(1) = -1;
rho_m    = 0.6;
for t = 2:T
    m_seq(t) = rho_m * m_seq(t-1);
end

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

% get direct PE impact on consumption

c_seq_PE = (C_ib + C_tau * dtaue_dib) * dib_dm * m_seq;

% get direct PE impact on labor supply

cs_seq_PE = (Cs_ib + Cs_tau * dtaue_dib) * dib_dm * m_seq;
l_seq_PE  = varphi/kappa_w * ((1+kappa_w) * zeros(T,1) - zeros(T,1) + zeros(T,1) ...
    - 1/R_b_SS * (zeros(T,1) - zeros(T,1) + zeros(T,1)) - kappa_w * gamma * cs_seq_PE);

%----------------------------------------------------------------
% Save Results
%----------------------------------------------------------------

scale = 1/max(abs(c_seq));

c_seq_mp      = scale*c_seq;
y_seq_mp      = scale*y_seq;
pi_seq_mp     = scale*pi_seq;
b_seq_mp      = scale*b_seq;
ib_seq_mp     = scale*ib_seq;
m_seq_mp      = scale*m_seq;
c_seq_mp_PE   = scale*c_seq_PE;
l_seq_mp_PE   = scale*l_seq_PE;

cd([path model task '/_results']);

save results_mp c_seq_mp y_seq_mp pi_seq_mp b_seq_mp ib_seq_mp rho_m m_seq_mp c_seq_mp_PE l_seq_mp_PE

cd([path model task]);