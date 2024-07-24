function C_y_pred = suffstats_fn(omega,theta,xi,psi,T);

% preparations

r_b_SS  = theta * (1/omega - 1 + xi)/(1/omega - 1) - 1;
omega_B = 1 - theta/(1 + r_b_SS);

% auxiliary BiU matrix

c0    = zeros(T,1);
c0(1) = omega_B;
for t = 2:T
    c0(t) = omega_B * theta^(t-1);
end

eta = zeros(T,1);
eta(1) = omega_B;
eta(2) = (psi*theta)/(1+r_b_SS) * eta(1);
for t = 3:T
    eta(t) = theta/(1+r_b_SS) * eta(t-1);
end

C_y_tilde = zeros(T,1);
C_y_tilde(:,1) = c0;
for t = 2:T
    C_y_tilde(1,t) = eta(t);
    C_y_tilde(2:T,t) = - eta(t) * (1+r_b_SS) * c0(1:T-1);
end

C_y_pred = zeros(T,T);
C_y_pred(:,1) = C_y_tilde(:,1);
for t = 2:T
    C_y_pred(:,t) = C_y_tilde(:,t) + [0;C_y_pred(1:T-1,t-1)];
end

% construct full matrix

mu = (omega - omega_B)/(1-omega_B);

C_y_pred = mu * eye(T) + (1-mu) * C_y_pred;