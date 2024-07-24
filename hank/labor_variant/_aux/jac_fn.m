function [M_c,M_cs,M_b] = jac_fn(Y_c,Y_cs,D);

%% PREPARATIONS: GLOBAL VARIABLES

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

% settings

global T step

% auxiliary variables

global c_opt_vec_SS Qt_big_SS

%% MAIN COMPUTATIONS

%----------------------------------------------------------------
% Construct P
%----------------------------------------------------------------

P_c = NaN(n_a*n_y,T-1);

for i_yT = 1:n_yT
    P_c(1+(i_yT-1)*(n_a*n_yP):i_yT*(n_a*n_yP),1) = c_opt_vec_SS(:,i_yT);
end

for t = 2:T-1
    P_c(:,t) = Qt_big_SS * P_c(:,t-1);
end

P_cs = NaN(n_a*n_y,T-1);
for i_yT = 1:n_yT
    P_cs(1+(i_yT-1)*(n_a*n_yP):i_yT*(n_a*n_yP),1) = (states(:,2) * grid_yT(i_yT)') .* c_opt_vec_SS(:,i_yT).^(-gamma);
end

for t = 2:T-1
    P_cs(:,t) = Qt_big_SS * P_cs(:,t-1);
end

P_b = NaN(n_a*n_y,T-1);
P_b(:,1) = repmat(grid_a',n_y,1);
for t = 2:T
    P_b(:,t) = Qt_big_SS * P_b(:,t-1);
end

%----------------------------------------------------------------
% Construct F
%----------------------------------------------------------------

F_c  = [Y_c;P_c'*D];
F_cs = [Y_cs;P_cs'*D];
F_b  = P_b'*D;

%----------------------------------------------------------------
% Construct M
%----------------------------------------------------------------

M_c = F_c;
for t = 2:T
    M_c(2:end,t) = M_c(2:end,t) + M_c(1:end-1,t-1);
end
M_c = M_c/C_SS;

M_cs = F_cs;
for t = 2:T
    M_cs(2:end,t) = M_cs(2:end,t) + M_cs(1:end-1,t-1);
end
M_cs = -1/gamma * M_cs/C_star_SS^(-gamma); % adjust to move derivative from C_star^(-gamma) to C_star

M_b = F_b;
for t = 2:T
    M_b(2:end,t) = M_b(2:end,t) + M_b(1:end-1,t-1);
end
M_b = M_b/B_SS;