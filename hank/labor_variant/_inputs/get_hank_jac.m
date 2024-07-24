%% HANK JACOBIANS
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
model = '/hank/labor_variant';
task = '/_inputs';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path model '/_aux']);
addpath([path model task '/_income_process']);
addpath([path model task '/_results']);
cd([path model task]);

%% PREPARATIONS

%----------------------------------------------------------------
% Global Variables
%----------------------------------------------------------------

% aggregate parameters

global beta beta_hat gamma probdeath varphi epsilon_w phi_w kappa_w wealth_0_pos ...
     epsilon_p phi_p ...
     labtax TransY_ratio BY_ratio rho_tr phi_pi phi_y phi_dy rho_dr
 
% household parameters

global a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap borr_wedge

% steady-state quantities

global C_SS C_star_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS r_b_SS B_SS D_SS Z_SS ...
    C_sd_SS wealth_pctls C_pctl_SS liqwealth_indic_SS lambda_pctl_SS n_pctls ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

% other quantities

global grid_a spliorder states states_div states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

load param_agg
load param_households
load SS
load aux

%----------------------------------------------------------------
% Time Horizon
%----------------------------------------------------------------

global T

T = 500;

%% GET ALL M MATRICES

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

global step

step = 5 * 10^(-3);

%----------------------------------------------------------------
% Prep Run
%----------------------------------------------------------------

global c_opt_vec_SS Qt_big_SS

w_seq     = zeros(T,1);
l_seq     = zeros(T,1);
trans_seq = zeros(T,1);
ib_seq    = zeros(T,1);
pi_seq    = zeros(T,1);
d_seq     = zeros(T,1);
zeta_seq  = zeros(T,1);
    
solve_hh_problem_lab

c_opt_vec_SS = c_opt_t(:,:,end);

Qt_SS     = Qt;
Qt_big_SS = NaN(n_y*n_a,n_y*n_a);

for i_yT = 1:n_yT
    Qt_big_SS(:,1+(i_yT-1)*(n_a*n_yP):i_yT*(n_a*n_yP)) = repmat(Qt_SS,n_yT,1) * yT_dist(i_yT);
end

%----------------------------------------------------------------
% Income
%----------------------------------------------------------------

vars = zeros(7,1);
vars(1) = step;

[Y_c_w,Y_cs_w,D_w] = YD_fn(vars);

[C_w,Cs_w,B_w] = jac_fn(Y_c_w,Y_cs_w,D_w);

C_l  = C_w;
Cs_l = Cs_w;
B_l  = B_w;

%----------------------------------------------------------------
% Transfers
%----------------------------------------------------------------

vars = zeros(7,1);
vars(3) = step;

[Y_c_tau,Y_cs_tau,D_tau] = YD_fn(vars);

[C_tau,Cs_tau,B_tau] = jac_fn(Y_c_tau,Y_cs_tau,D_tau);

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

vars = zeros(7,1);
vars(5) = step;

[Y_c_pi,Y_cs_pi,D_pi] = YD_fn(vars);

[C_pi,Cs_pi,B_pi] = jac_fn(Y_c_pi,Y_cs_pi,D_pi);

%----------------------------------------------------------------
% Nominal Interest Rate
%----------------------------------------------------------------

C_ib   = [-C_pi(:,2:T),[0;-C_pi(1:T-1,T)]];
Cs_ib  = [-Cs_pi(:,2:T),[0;-Cs_pi(1:T-1,T)]];
B_ib   = [-B_pi(:,2:T),[0;-B_pi(1:T-1,T)]];

%----------------------------------------------------------------
% Dividend
%----------------------------------------------------------------

vars = zeros(7,1);
vars(6) = step;

[Y_c_d,Y_cs_d,D_d] = YD_fn(vars);

[C_d,Cs_d,B_d] = jac_fn(Y_c_d,Y_cs_d,D_d);

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

vars = zeros(7,1);
vars(7) = step;

[Y_c_b,Y_cs_b,D_b] = YD_fn(vars);

[C_b,Cs_b,B_b] = jac_fn(Y_c_b,Y_cs_b,D_b);

%% SAVE RESULTS

cd([path model task '/_results']);

save jac_matrices C_w Cs_w B_w C_l Cs_l B_l C_ib Cs_ib B_ib C_pi Cs_pi B_pi ...
    C_tau Cs_tau B_tau C_d Cs_d B_d C_b Cs_b B_b step

cd([path model task]);