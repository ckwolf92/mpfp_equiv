%% HANK STEADY STATE
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
model = '/hank/baseline';
task = '/_inputs';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path model '/_aux']);
addpath([path model task '/_income_process']);
cd([path model task]);

%% CALIBRATION TARGETS

% set wealth calibration

lw_base_indic = 1;
lw_low_indic  = 0;
lw_high_indic = 0;

% list of targets

AY_target             = (2.92 - 0.26) * 4;
if lw_base_indic == 1
    BY_target = 1.5;
elseif lw_low_indic == 1
    BY_target = 1.5/3;
elseif lw_high_indic == 1
    BY_target = 1.5 * 5;
end
TransY_target         = 0.05;
tax_target            = 0.3;
profshare_target      = 0.06;
r_b_target            = 0.01;

%% ECONOMIC PARAMETERS

%----------------------------------------------------------------
% Household
%----------------------------------------------------------------

% preferences

beta      = [];
gamma     = 1;
probdeath = 1/180;

varphi    = 0.5;
epsilon_w = 10;
lambda_w  = 1/(epsilon_w - 1);
phi_w     = 0.7;

% wealth lower bound

a_lb = 0;

% income risk

load logyPgrid.txt
load yPtrans.txt
grid_yP = exp(logyPgrid);
n_yP = length(grid_yP);
Pi_yP = yPtrans;
clear logyPgrid yPtrans

load logyTgrid.txt
load yTdist.txt
grid_yT = exp(logyTgrid);
n_yT = length(grid_yT);
Pi_yT = repmat(yTdist',n_yT,1);
clear logyTgrid yTdist

n_y = n_yT * n_yP;
grid_y = repmat(grid_yP,1,n_yT) .* reshape(repmat(grid_yT,1,n_yP)',1,n_y);
indx_yP = reshape(repmat((1:1:n_yP)',1,n_yT),n_y,1)';
indx_yT = reshape(repmat(1:1:n_yT,n_yP,1),n_y,1)';

Pi_y = NaN(n_y,n_y);
for i_y = 1:n_y
    for i_yy = 1:n_y
        Pi_y(i_y,i_yy) = Pi_yP(indx_yP(i_y),indx_yP(i_yy)) * Pi_yT(indx_yT(i_y),indx_yT(i_yy));
    end
end

y_dist  = ergodicdist(Pi_y)'; 
yT_dist = ergodicdist(Pi_yT)';
yP_dist = ergodicdist(Pi_yP)';

mean_y_raw = sum(grid_y(:).*y_dist(:));
grid_y     = grid_y/mean_y_raw;
grid_yP    = grid_yP/mean_y_raw;
mean_y     = grid_y * y_dist';

% dividends

div_het = 1;

if div_het == 1

    chi = 3;

    yP_median_indx = find(cumsum(yP_dist)>0.5,1);
    div_share      = 0 * grid_yP;
    for i_yP = yP_median_indx+1:n_yP
        div_share(i_yP) = (grid_yP(i_yP) - grid_yP(yP_median_indx))^chi;
    end
    div_share = div_share./(div_share * yP_dist');

else
    
    div_share = ones(1,n_yP);

end

%----------------------------------------------------------------
% Production
%----------------------------------------------------------------

epsilon_p  = 1/profshare_target;
phi_p      = 0.85;
AY_ratio   = AY_target;

lambda_p   = 1/(epsilon_p - 1);

%----------------------------------------------------------------
% Government
%----------------------------------------------------------------

labtax       = tax_target;
BY_ratio     = BY_target;
TransY_ratio = TransY_target;

rho_tr       = 0;
phi_pi       = 1.5;
phi_y        = 0;
phi_dy       = 0;

rho_dr       = 0.85;

%% STEADY-STATE CALIBRATION

% interest rate

r_b_SS         = r_b_target;
annuity_gap    = probdeath/(1-probdeath);
borr_wedge     = 0;

% total output

Z_SS = 1;
Y_SS = epsilon_p/(epsilon_p-1);

% marginal cost for retailers

P_I_SS = (epsilon_p-1)/epsilon_p;

% total labor income

L_SS = Y_SS/Z_SS;
W_SS = P_I_SS * Y_SS/L_SS;

% total debt

B_SS = BY_ratio * Y_SS;

% total profits

D_SS = Y_SS - W_SS * L_SS;

% value of the corporate sector

V_SS = D_SS/(1-1/(1 + r_b_SS));
A_SS = V_SS - D_SS;

% total government transfers and spending

Trans_SS = TransY_ratio * Y_SS;
G_SS     = labtax * W_SS * L_SS + B_SS - Trans_SS - (1+r_b_SS) * B_SS;

% transfer for MPC distribution

gift = 10^(-3) * Trans_SS;

%% SOLUTION PARAMETERS

%----------------------------------------------------------------
% Grids
%----------------------------------------------------------------

n_a        = 80; % number of grid points 
a_min      = a_lb;
a_max      = 30 * max(grid_yP) * W_SS * L_SS;
spliorder  = [1 1];

%----------------------------------------------------------------
% Beta Loop
%----------------------------------------------------------------

beta_init   = 0.99 * 1/(1 + r_b_SS);
beta_guess  = beta_init;
beta_upd    = 0.4;
beta_tol    = 10^(-4);
beta_it_max = 50;
beta_ub     = 1;
beta_lb     = 0.90;

%----------------------------------------------------------------
% EGP Iteration
%----------------------------------------------------------------

EGP_tol      = 10^(-8);
disp_EGPdist = 0;

%% ASSEMBLE GRID

%----------------------------------------------------------------
% Asset Grid
%----------------------------------------------------------------

if a_min >= 0
    
    gridparam_a = 0.3; %linear = 1, 1 L-shaped = 0;

    grid_a = linspace(0,1,n_a);
    grid_a = grid_a.^(1./gridparam_a);
    grid_a = a_min + (a_max-a_min).*grid_a;
    
else

    n_a_pos = round(3/4 * n_a);
    n_a_neg = n_a - n_a_pos;

    gridparam_a_pos = 0.3; %linear = 1, 1 L-shaped = 0;
    grid_a_pos = linspace(0,1,n_a_pos);
    grid_a_pos = grid_a_pos.^(1./gridparam_a_pos);
    grid_a_pos = 0 + (a_max-0).*grid_a_pos;

    gridparam_a_neg = 0.7; %linear = 1, 1 L-shaped = 0;
    grid_a_neg = linspace(0,1,n_a_neg+1);
    grid_a_neg = grid_a_neg.^(1./gridparam_a_neg);
    grid_a_neg = sort(-(0 + (-a_min-0).*grid_a_neg),'ascend');

    grid_a = [grid_a_neg(1:end-1),grid_a_pos];
    
end

wealth_0_pos = find(grid_a == 0);

%----------------------------------------------------------------
% Splines
%----------------------------------------------------------------

% put productivity and asset grids together

n_s  = n_a * n_yP;

fspace          = fundef({'spli',grid_a,0,spliorder(1)},...
                         {'spli',grid_yP,0,spliorder(2)});
states_grid  = funnode(fspace);
states       = gridmake(states_grid);

states_a   = states(:,1);
states_yP  = states(:,2);

Phi_yP     = splibas(grid_yP,0,spliorder(2),states(:,2));
Phi_A      = splibas(grid_a,0,spliorder(1),states(:,1));
Phi        = dprod(Phi_yP,Phi_A);
Emat_yP    = kron(Pi_yP,speye(n_a));

%----------------------------------------------------------------
% Return Grid
%----------------------------------------------------------------

r_b_grid = (states_a >= 0) .* (r_b_SS + annuity_gap) + (states_a < 0) .* (r_b_SS + annuity_gap + borr_wedge);

%----------------------------------------------------------------
% Dividend Grid
%----------------------------------------------------------------

states_div   = reshape(repmat(div_share,n_a,1),n_yP*n_a,1);

%% MAIN SOLUTION LOOP

for beta_it = 1:beta_it_max
    
%----------------------------------------------------------------
% Discount Facor
%----------------------------------------------------------------

beta     = beta_guess;
beta_hat = beta * (1 - probdeath);  
 
%----------------------------------------------------------------
% Endogenous Gridpoint Iteration
%----------------------------------------------------------------

% preparations

dist_EGP = 1;
EGP_it   = 0;

% initial guess

cp_opt  = (1-labtax) * W_SS * L_SS * (states_yP * grid_yT') + r_b_grid .* states_a + Trans_SS + states_div * D_SS;
mutilde_opt = Emat_yP * (((1 + r_b_grid) .* cp_opt.^(-gamma)) * yT_dist');

% iteration

while dist_EGP > EGP_tol
    
% one step for EGP

[c_opt,ap_opt,mutilde_upd] = EGP_fun(mutilde_opt,r_b_grid,W_SS,L_SS,Trans_SS + states_div * D_SS,0,...
    beta_hat,gamma,labtax,states,grid_a,Emat_yP,grid_yT,yT_dist);
    
% update

dist_EGP = norm(mutilde_upd - mutilde_opt)/norm(mutilde_opt);

if disp_EGPdist == 1 && mod(EGP_it,100) == 0
    disp(dist_EGP)
end

mutilde_opt = mutilde_upd;
EGP_it      = EGP_it + 1;

end

mutilde_SS = mutilde_opt;

%----------------------------------------------------------------
% Distribution
%----------------------------------------------------------------

ap_opt     = max(min(ap_opt,a_max),a_min);
fspaceerga = fundef({'spli',grid_a,0,1});

QZ_live = kron(Pi_yP,ones(n_a,1));
QA_live = 0;
for i_yT = 1:n_yT
    QA_live = QA_live + yT_dist(i_yT) * funbas(fspaceerga,ap_opt(:,i_yT));
end
Q_live  = dprod(QZ_live,QA_live);

QA_death = sparse(n_s,n_a);
QA_death(:,wealth_0_pos) = 1;
QZ_death = repmat(yP_dist,n_s,1);
Q_death  = dprod(QZ_death,QA_death);

Q = (1-probdeath) * Q_live + probdeath * Q_death;

lambda_SS      = ergodicdist(Q,2);
lambda_vec_SS  = lambda_SS;
lambda_SS      = permute(reshape(lambda_SS,[n_a,n_yP]),[2 1]);
lambdafull_SS  = kron(yT_dist',lambda_SS);

%----------------------------------------------------------------
% Compute Aggregates
%----------------------------------------------------------------

% re-order policy function

c_opt_SS  = permute(reshape(c_opt,[n_a,n_y]),[2 1]);
ap_opt_SS = permute(reshape(ap_opt,[n_a,n_y]),[2 1]);

% compute other aggregates

C_grid = c_opt_SS;
C_SS   = sum(sum(C_grid(:,:).*lambdafull_SS(:,:)));

A_grid    = repmat(grid_a,n_y,1);
Omega_SS  = sum(sum(A_grid(:,:).*lambdafull_SS(:,:)));

%----------------------------------------------------------------
% Check Asset Market Clearing
%----------------------------------------------------------------

Wealth_err = Omega_SS - B_SS;

% ----------------------------------------------------------------
% Beta Update
% ----------------------------------------------------------------

if Wealth_err < -beta_tol
    disp(['Discount factor: ' num2str(beta_guess) ', Discount factor too low: ' num2str(Wealth_err) ]);
    beta_lb = beta_guess;
    beta_guess = (1-beta_upd) * beta_guess + beta_upd * beta_ub;
elseif Wealth_err > beta_tol
    disp(['Discount factor: ' num2str(beta_guess) ', Discount factor too high: ' num2str(Wealth_err) ]);
    beta_ub = beta_guess;
    beta_guess = (1-beta_upd) * beta_guess + beta_upd * beta_lb;
elseif abs(Wealth_err) <= beta_tol
    disp(['Steady State Found, Discount factor = ' num2str(beta_guess)]);
    break
end

end

%% POST-COMPUTATIONS

% other aggregates

Pi_SS  = 1;
R_b_SS = 1 + r_b_SS;
R_n_SS = R_b_SS;

% NKPC slopes

kappa_w   = (((1-phi_w*1/R_b_SS)*(1-phi_w))/(phi_w*(1+1/varphi*(1+1/lambda_w))));
kappa_p   = 1/((1-1/R_b_SS*phi_p)*(1-phi_p)/(phi_p));

% consumption dispersion

C_sd_SS = sqrt(sum(sum((c_opt - C_SS).^2 .* kron(lambda_vec_SS,yT_dist))));

%% STEADY-STATE STATISTICS: COMPUTATIONS

%----------------------------------------------------------------
% MPC
%----------------------------------------------------------------

[c_opt_MPC,~,~] = EGP_fun(mutilde_opt,r_b_grid,W_SS,L_SS,Trans_SS + states_div * D_SS + gift,0,...
    beta_hat,gamma,labtax,states,grid_a,Emat_yP,grid_yT,yT_dist);

c_opt_MPC = permute(reshape(c_opt_MPC,[n_a,n_y]),[2 1]);
MPC_dist_SS = (c_opt_MPC - c_opt_SS)./gift;
MPC_avg_SS = sum(MPC_dist_SS(:) .* lambdafull_SS(:));

%----------------------------------------------------------------
% Distributional Aggregates
%----------------------------------------------------------------

liqwealth_dist_SS = sum(lambda_SS,1);

wealth_pctls = linspace(0,1,15);
n_pctls      = length(wealth_pctls)-1;

C_pctl_SS          = NaN(1,n_pctls);
liqwealth_indic_SS = NaN(size(lambdafull_SS,1),size(lambdafull_SS,2),n_pctls);
lambda_pctl_SS     = NaN(size(lambdafull_SS,1),size(lambdafull_SS,2),n_pctls);

for i_pctl = 1:n_pctls
    [C_pctl_SS(1,i_pctl),liqwealth_indic_SS(:,:,i_pctl),lambda_pctl_SS(:,:,i_pctl)] = wealth_pctl_fn(wealth_pctls(i_pctl),wealth_pctls(i_pctl+1),...
                                                        liqwealth_dist_SS,lambdafull_SS,C_grid,grid_a);
end

%% SAVE RESULTS

cd([path model task '/_results']);

if lw_base_indic == 1

% aggregate parameters

save param_agg beta beta_hat gamma probdeath varphi epsilon_w phi_w kappa_w wealth_0_pos ...
     epsilon_p phi_p kappa_p ...
     labtax TransY_ratio BY_ratio rho_tr phi_pi phi_y phi_dy rho_dr

% household parameters

save param_households a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap borr_wedge

% steady state

save SS C_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS r_b_SS B_SS D_SS Z_SS ...
    C_sd_SS wealth_pctls C_pctl_SS liqwealth_indic_SS lambda_pctl_SS n_pctls ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

% other quantities

save aux grid_a spliorder states states_div states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

elseif lw_low_indic == 1

% aggregate parameters

save param_agg_low beta beta_hat gamma probdeath varphi epsilon_w phi_w kappa_w wealth_0_pos ...
     epsilon_p phi_p kappa_p ...
     labtax TransY_ratio BY_ratio rho_tr phi_pi phi_y phi_dy rho_dr

% household parameters

save param_households_low a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap borr_wedge

% steady state

save SS_low C_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS r_b_SS B_SS D_SS Z_SS ...
    C_sd_SS wealth_pctls C_pctl_SS liqwealth_indic_SS lambda_pctl_SS n_pctls ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

% other quantities

save aux_low grid_a spliorder states states_div states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

elseif lw_high_indic == 1

% aggregate parameters

save param_agg_high beta beta_hat gamma probdeath varphi epsilon_w phi_w kappa_w wealth_0_pos ...
     epsilon_p phi_p kappa_p ...
     labtax TransY_ratio BY_ratio rho_tr phi_pi phi_y phi_dy rho_dr

% household parameters

save param_households_high a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap borr_wedge

% steady state

save SS_high C_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS r_b_SS B_SS D_SS Z_SS ...
    C_sd_SS wealth_pctls C_pctl_SS liqwealth_indic_SS lambda_pctl_SS n_pctls ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

% other quantities

save aux_high grid_a spliorder states states_div states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

end

cd([path model task]);