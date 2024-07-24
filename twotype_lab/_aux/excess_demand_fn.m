function excess_demand = excess_demand_fn(guess_seq);

%% GLOBAL VARIABLES

global varphi delta chi_B chi_H epsilon_p phi_p kappa_p ...
     labtax rho_tr phi_pi phi_y phi_dy rho_dr

global beta gamma

global C_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS B_SS D_SS Z_SS r_b_SS

global C_w B_w L_w C_tau B_tau L_tau C_ib B_ib L_ib C_pi B_pi L_pi C_d B_d L_d C_b B_b L_b

global T

global dib_dpi dib_dy dib_dm ...
    dtaue_dtaux dbg_dtaux dtaue_dw dbg_dw dtaue_dl dbg_dl dtaue_dpi dbg_dpi dtaue_dib dbg_dib

global taux_seq m_seq zeta_seq

%% COLLECT INPUTS

pi_seq  = guess_seq(1:T,1);
l_seq   = guess_seq(T+1:2*T,1);

%% GET OUTCOMES

get_aggregates

%% CHECK ACCURACY

excess_demand_1 = b_seq - bg_seq;
excess_demand_2 = l_seq - l_seq_supply;

excess_demand = [excess_demand_1;excess_demand_2];