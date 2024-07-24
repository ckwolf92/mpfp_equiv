function [c,b,l] = c_B_fn(exo,A_inv);

global labtax W_SS L_SS Trans_B_SS r_b_SS B_B_SS D_SS gamma beta

% get inputs

w_hat     = exo.w_hat;
tau_hat   = exo.tau_hat;
ib_hat    = exo.ib_hat;
pi_hat    = exo.pi_hat;
d_hat     = exo.d_hat;
zeta_hat  = exo.zeta_hat;

T = length(w_hat);

% solve equilibrium system

b = [(1-labtax) * W_SS * L_SS * w_hat + Trans_B_SS * tau_hat + (1 + r_b_SS) * B_B_SS * ([0;ib_hat(1:end-1)] - pi_hat) + D_SS * d_hat; ...
        -1/gamma * beta * (1+r_b_SS) * (ib_hat - [pi_hat(2:end);0]) + zeta_hat; ...
        (1-labtax) * W_SS * w_hat; ...
        zeros(T,1)];

sol = A_inv * b;

c = sol(1:T);
b = sol(T+1:2*T);
l = sol(2*T+1:3*T);

end