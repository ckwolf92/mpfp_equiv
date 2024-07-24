%% GET OUTCOMES

% output

y_seq = l_seq;

% monetary block

ib_seq = dib_dpi * pi_seq + dib_dy * y_seq + dib_dm * m_seq;

% intermediate goods price

p_I_seq = kappa_p * (pi_seq - 1/(1+r_b_SS) * [pi_seq(2:end);0]);

% wage

w_seq = p_I_seq;

% dividends

d_seq = 1/D_SS * (Y_SS * y_seq - W_SS * L_SS * (w_seq + l_seq));

% transfers and debt

bg_seq   = dbg_dtaux * taux_seq + dbg_dw * w_seq + dbg_dl * l_seq + dbg_dpi * pi_seq + dbg_dib * ib_seq;
taue_seq = dtaue_dtaux * taux_seq + dtaue_dw * w_seq + dtaue_dl * l_seq + dtaue_dpi * pi_seq + dtaue_dib * ib_seq;
tau_seq  = taux_seq + taue_seq;

% consumption

c_seq = C_w * w_seq + C_tau * tau_seq + C_ib * ib_seq + C_pi * pi_seq + C_d * d_seq + C_b * zeta_seq;
b_seq = B_w * w_seq + B_tau * tau_seq + B_ib * ib_seq + B_pi * pi_seq + B_d * d_seq + B_b * zeta_seq;

% labor supply

l_seq_supply = L_w * w_seq + L_tau * tau_seq + L_ib * ib_seq + L_pi * pi_seq + L_d * d_seq + L_b * zeta_seq;