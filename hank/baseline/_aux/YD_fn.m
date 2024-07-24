function [Y_c,D] = YD_fn(vars);

%% PREPARATIONS: GLOBAL VARIABLES

% aggregate parameters

global beta beta_hat gamma probdeath varphi epsilon_w phi_w kappa_w wealth_0_pos ...
     epsilon_p phi_p kappa_p ...
     labtax TransY_ratio BY_ratio rho_tr phi_pi phi_y phi_dy rho_dr
 
% household parameters

global a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap borr_wedge

% steady-state quantities

global C_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS r_b_SS B_SS D_SS Z_SS ...
    C_sd_SS wealth_pctls C_pctl_SS liqwealth_indic_SS lambda_pctl_SS n_pctls ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

% other quantities

global grid_a spliorder states states_div states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

% settings

global T step

%% MAIN COMPUTATIONS

%----------------------------------------------------------------
% Price Sequences
%----------------------------------------------------------------

% input sequences

w_seq     = zeros(T,1);
l_seq     = zeros(T,1);
trans_seq = zeros(T,1);
ib_seq    = zeros(T,1);
pi_seq    = zeros(T,1);
d_seq     = zeros(T,1);
zeta_seq  = zeros(T,1);

w_seq(end)     = vars(1);
l_seq(end)     = vars(2);
trans_seq(end) = vars(3);
ib_seq(end-1)  = vars(4);
pi_seq(end)    = vars(5);
d_seq(end)     = vars(6);
zeta_seq(end)  = vars(7);

% get back to levels

r_b_grid_seq  = NaN(n_s,T);
r_b_grid_seq(:,1) = (states_a >= 0) .* (R_b_SS * exp(0-pi_seq(1)) - 1 + annuity_gap) ...
    + (states_a < 0) .* (R_b_SS * exp(0-pi_seq(1)) - 1 + annuity_gap + borr_wedge);
for t = 2:T
    r_b_grid_seq(:,t) = (states_a >= 0) .* (R_b_SS * exp(ib_seq(t-1)-pi_seq(t)) - 1 + annuity_gap) ...
        + (states_a < 0) .* (R_b_SS * exp(ib_seq(t-1)-pi_seq(t)) - 1 + annuity_gap + borr_wedge);
end
W_seq         = W_SS * exp(w_seq);
L_seq         = L_SS * exp(l_seq);
Trans_seq     = Trans_SS * exp(trans_seq);
D_seq         = D_SS * exp(d_seq);

%----------------------------------------------------------------
% Placeholders
%----------------------------------------------------------------

c_opt_t   = zeros(n_s,n_yT,T);
ap_opt_t  = zeros(n_s,n_yT,T);

mutilde_t        = zeros(n_s,T+1);
mutilde_t(:,T+1) = mutilde_SS;

%----------------------------------------------------------------
% Get Policies
%----------------------------------------------------------------

for t = T:-1:1
    [c_opt_t(:,:,t),ap_opt_t(:,:,t),mutilde_t(:,t)] = EGP_fun(mutilde_t(:,t+1),r_b_grid_seq(:,t),W_seq(t),L_seq(t),Trans_seq(t) + states_div * D_seq(t),zeta_seq(t),...
                            beta_hat,gamma,labtax,states,grid_a,Emat_yP,grid_yT,yT_dist);
end

%----------------------------------------------------------------
% Get Distribution
%----------------------------------------------------------------

lambda_seq = NaN(n_yP,n_a,T+1);
lambdafull_seq = NaN(n_y,n_a);
QZ_live = kron(Pi_yP,ones(n_a,1));

QAt_death = sparse(n_s,n_a);
QAt_death(:,wealth_0_pos) = 1;
QZt_death = repmat(yP_dist,n_s,1);
Qt_death  = dprod(QZt_death,QAt_death);

Q_list = cell(T,1);

for t = 1:T+1
    if t == 1
        lambda_aux = lambda_vec_SS(:);
        lambda_lag = lambda_aux;
        lambda_aux = permute(reshape(lambda_aux,[n_a,n_yP]),[2 1]);  
        lambda_seq(:,:,t) = lambda_aux;
    else
        ap_opt     = max(min(ap_opt_t(:,:,t-1),a_max),a_min);
        fspaceerga = fundef({'spli',grid_a,0,1});
        QAt_live = 0;
        for i_yT = 1:n_yT
            QAt_live = QAt_live + yT_dist(i_yT) * funbas(fspaceerga,ap_opt(:,i_yT));
        end
        Qt_live  = dprod(QZ_live,QAt_live);
        
        Qt = (1-probdeath) * Qt_live + probdeath * Qt_death;
        
        % update distribution
        
        lambda_aux = Qt' * lambda_lag;
        lambda_lag = lambda_aux;
        lambda_aux = permute(reshape(lambda_aux,[n_a,n_yP]),[2 1]);
        lambda_seq(:,:,t) = lambda_aux;
        Q_list{t-1} = Qt;
        
    end
    
    lambdafull_seq(:,:,t) = kron(yT_dist',lambda_seq(:,:,t));
    
end

%----------------------------------------------------------------
% Get Y & D
%----------------------------------------------------------------

Y_c = NaN(1,T);
D = NaN(n_a*n_y,T);

for t = 1:T
    
c_tt = permute(reshape(c_opt_t(:,:,T-t+1),[n_a,n_y]),[2 1]);

Y_c(t) = (sum(sum(c_tt .* lambdafull_SS)) - C_SS)/step;

for i_yT=1:n_yT

    D(1+(i_yT-1)*(n_a*n_yP):i_yT*(n_a*n_yP),t) = yT_dist(i_yT) * (Q_list{T-t+1}' * lambda_vec_SS - lambda_vec_SS)/step;

end

end