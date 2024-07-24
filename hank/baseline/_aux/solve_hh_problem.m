%----------------------------------------------------------------
% Price Sequences
%----------------------------------------------------------------

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
        
    end
    
    lambdafull_seq(:,:,t) = kron(yT_dist',lambda_seq(:,:,t));
    
end

%----------------------------------------------------------------
% Aggregation
%----------------------------------------------------------------

% household-level output sequences

B_grid_t        = NaN(n_y,n_a,T+1);
C_grid_t        = NaN(n_y,n_a,T+1);
for t = 1:T
    B_grid_t(:,:,t) = repmat(grid_a,n_y,1);
    C_grid_t(:,:,t) = permute(reshape(c_opt_t(:,:,t),[n_a,n_y]),[2 1]);
end
B_grid_t(:,:,T+1) = repmat(grid_a,n_y,1);
C_grid_t(:,:,T+1) = permute(reshape(c_opt_t(:,:,T),[n_a,n_y]),[2 1]);

c_opt_t(:,:,T+1) = c_opt_t(:,:,T);

% consumption & assets

B_seq_HH     = NaN(T+1,1);
b_seq_HH     = NaN(T+1,1);
C_seq        = NaN(T+1,1);
c_seq        = NaN(T+1,1);
C_pctl_seq   = NaN(T+1,n_pctls);
c_pctl_seq   = NaN(T+1,n_pctls);
C_sd_seq     = NaN(T+1,1);
c_sd_seq     = NaN(T+1,1);

for t = 1:T+1
    B_seq_HH(t)     = sum(sum(B_grid_t(:,:,t).*lambdafull_seq(:,:,t)));   
    C_seq(t)        = sum(sum(C_grid_t(:,:,t).*lambdafull_seq(:,:,t)));
    C_sd_seq(t)     = sqrt(sum(sum((C_grid_t(:,:,t) - C_seq(t)).^2 .* lambdafull_seq(:,:,t))));
    
    for i_pctl = 1:n_pctls
        lambdafull_pctl_seq_t = liqwealth_indic_SS(:,:,i_pctl) .* lambdafull_SS;
        lambdafull_pctl_seq_t = lambdafull_pctl_seq_t ./ sum(lambdafull_pctl_seq_t(:));
        C_pctl_seq(t,i_pctl) = sum(sum(C_grid_t(:,:,t).*lambdafull_pctl_seq_t));
    end
    
    b_seq_HH(t)     = log(B_seq_HH(t)/B_SS);
    c_seq(t)        = log(C_seq(t)/C_SS);
    c_sd_seq(t)     = log(C_sd_seq(t)/C_sd_SS);
    for i_pctl = 1:n_pctls
        c_pctl_seq(t,i_pctl) = log(C_pctl_seq(t,i_pctl)/C_pctl_SS(i_pctl));
    end
end

c_seq      = c_seq(1:T,1);
c_pctl_seq = c_pctl_seq(1:T,:);
b_seq_HH   = b_seq_HH(2:T+1,1);