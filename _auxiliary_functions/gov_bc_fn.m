function [tau_e_seq,b_g_seq] = gov_bc_fn(taux_seq,ib_seq,pi_seq,w_seq,l_seq);

global T rho_dr r_b_SS labtax W_SS L_SS B_SS Trans_SS

b_g_seq = zeros(T,1);

b_g_seq(1) = 1/B_SS * (Trans_SS * taux_seq(1) + (1 + r_b_SS) * B_SS * (0 - pi_seq(1)) - labtax * W_SS * L_SS * (w_seq(1) + l_seq(1)));
for t = 2:T
    b_g_seq(t) = rho_dr * b_g_seq(t-1) + 1/B_SS * (Trans_SS * taux_seq(t) + ...
        (1 + r_b_SS) * B_SS * (ib_seq(t-1) - pi_seq(t)) - labtax * W_SS * L_SS * (w_seq(t) + l_seq(t)));
end

tau_seq = 1/Trans_SS * (labtax * W_SS * L_SS * (w_seq + l_seq) + B_SS * b_g_seq - (1 + r_b_SS) * B_SS * [0;ib_seq(1:end-1)] ...
    + (1 + r_b_SS) * B_SS * pi_seq - (1 + r_b_SS) * B_SS * [0;b_g_seq(1:end-1)]);
tau_e_seq = tau_seq - taux_seq;