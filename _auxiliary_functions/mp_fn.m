function ib_seq = mp_fn(pi_seq,y_seq,m_seq);

global T rho_tr phi_pi phi_y phi_dy

ib_seq = zeros(T,1);

ib_seq(1) = (1-rho_tr) * phi_pi * pi_seq(1) + (1-rho_tr) * (phi_y + phi_dy) * y_seq(1) - (1-rho_tr) * phi_dy * 0 + m_seq(1);
for t = 2:T
    ib_seq(t) = rho_tr * ib_seq(t-1) + (1-rho_tr) * phi_pi * pi_seq(t) + (1-rho_tr) * (phi_y + phi_dy) * y_seq(t) ...
        - (1-rho_tr) * phi_dy * y_seq(t-1) + m_seq(t);
end