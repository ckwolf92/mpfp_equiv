%% GENERATE INPUTS FOR TWO-TYPE MPE ANALYSIS
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
model = '/twotype_lab';
task = '/_inputs';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path model '/_aux']);
cd([path model task]);

%% PARAMETERS

%----------------------------------------------------------------
% Preferences
%----------------------------------------------------------------

global beta gamma eta alpha delta_B delta_H mu varphi chi_B chi_H

gamma   = 1;
eta     = 0.5;
alpha   = 0.1;
delta_B = 0.03;
delta_H = 0.95;
varphi  = 0.5;

mu    = 0.25;

%----------------------------------------------------------------
% Firms
%----------------------------------------------------------------

epsilon_p  = 1/0.06;
phi_p      = 0.85;
lambda_p   = 1/(epsilon_p - 1);

%----------------------------------------------------------------
% Government Policy
%----------------------------------------------------------------

rho_tr       = 0;
phi_pi       = 1.5;
phi_y        = 0;
phi_dy       = 0;

rho_dr       = 0.85;

%----------------------------------------------------------------
% Steady-State Objects
%----------------------------------------------------------------

global labtax Z_SS Y_SS P_I_SS B_SS r_b_SS Trans_SS L_SS W_SS G_SS C_SS D_SS ...
    B_B_SS B_H_SS C_B_SS C_H_SS L_B_SS L_H_SS X_B_SS X_H_SS Trans_H_SS Trans_B_SS ...
    Pi_SS R_b_SS R_n_SS

labtax   = 0.3;
Z_SS     = 1;
Y_SS     = epsilon_p/(epsilon_p-1);
P_I_SS   = (epsilon_p-1)/epsilon_p;
B_SS     = 1.5 * Y_SS;
r_b_SS   = 0.01;
Trans_SS = 0.05 * Y_SS;
L_SS     = Y_SS/Z_SS;
W_SS     = P_I_SS * Y_SS/L_SS;
G_SS     = labtax * W_SS * L_SS + B_SS - Trans_SS - (1+r_b_SS) * B_SS;
C_SS     = Y_SS - G_SS;
D_SS     = Y_SS - W_SS * L_SS;

B_H_SS     = B_SS;
B_B_SS     = (B_SS - mu * B_H_SS)/(1-mu);
L_H_SS     = L_SS;
L_B_SS     = (L_SS - mu * L_H_SS)/(1-mu);
C_H_SS     = C_SS;
C_B_SS     = (C_SS - mu * C_H_SS)/(1-mu);
Trans_H_SS = C_H_SS - ((1-labtax) * W_SS * L_H_SS + D_SS + r_b_SS * B_H_SS);
Trans_B_SS = C_B_SS - ((1-labtax) * W_SS * L_B_SS + D_SS + r_b_SS * B_B_SS);

chi_B_fn_aux = @(chi) (chi * delta_B * (C_B_SS - chi * delta_B * (L_B_SS^(1 + 1/varphi))/(1+1/varphi))^(-gamma) ...
    + chi * (1-delta_B)) * L_B_SS^(1/varphi) ...
    - (1-labtax) * (C_B_SS - chi * delta_B * (L_B_SS^(1+1/varphi))/(1+1/varphi))^(-gamma) * W_SS;
chi_B = fsolve(chi_B_fn_aux,1);
X_B_SS = C_B_SS - chi_B * delta_B * (L_B_SS^(1 + 1/varphi))/(1+1/varphi);

chi_H_fn_aux = @(chi) (chi * delta_H * (C_H_SS - chi * delta_H * (L_H_SS^(1 + 1/varphi))/(1+1/varphi))^(-gamma) ...
    + chi * (1-delta_H)) * L_H_SS^(1/varphi) ...
    - (1-labtax) * (C_H_SS - chi * delta_H * (L_H_SS^(1+1/varphi))/(1+1/varphi))^(-gamma) * W_SS;
chi_H = fsolve(chi_H_fn_aux,1);
X_H_SS = C_H_SS - chi_H * delta_H * (L_H_SS^(1 + 1/varphi))/(1+1/varphi);

beta_fn_aux = @(beta) X_B_SS^(-gamma) * (1 - beta * (1 + r_b_SS)) - alpha * beta * B_B_SS^(-eta);
beta = fsolve(beta_fn_aux,0);

Pi_SS  = 1;
R_b_SS = 1 + r_b_SS;
R_n_SS = R_b_SS;

kappa_p   = 1/((1-1/R_b_SS*phi_p)*(1-phi_p)/(phi_p));

%% SETTINGS

global T

T    = 500;
step = 1;

exo.w_hat     = zeros(T,1);
exo.tau_hat   = zeros(T,1);
exo.ib_hat    = zeros(T,1);
exo.pi_hat    = zeros(T,1);
exo.d_hat     = zeros(T,1);
exo.zeta_hat  = zeros(T,1);

r_NPV = zeros(T,1);
for t = 1:T
    r_NPV(t) = (1/(1+r_b_SS))^(t-1);
end

%% DERIVATIVE MATRIVES: B TYPES

%----------------------------------------------------------------
% Auxiliary Matrix
%----------------------------------------------------------------

A = zeros(4*T,4*T); % order: (c, b, l, x)
for t = 1:T % budget constraint
    A(t,t) = C_B_SS;
    A(t,T+t) = B_B_SS;
    if t > 1
        A(t,T+t-1) = -(1+r_b_SS) * B_B_SS;
    end
    A(t,2*T+t) = -(1-labtax) * W_SS * L_B_SS;
end
for t = T+1:2*T % Euler equation
    if t < 2*T
        A(t,t+2*T) = 1;
        A(t,t+2*T+1) = - beta * (1+r_b_SS);
        A(t,t) = -eta/gamma * (1 - beta * (1+r_b_SS));
    else
        A(t,t+2*T) = 1;
        A(t,t) = 1;
    end
end
for t = 2*T+1:3*T % labor supply relation
    A(t,t) = 1/varphi * chi_B * delta_B * L_B_SS^(1/varphi) + chi_B * (1-delta_B) * L_B_SS^(1/varphi) * X_B_SS^gamma * 1/varphi;
    A(t,t+T) = chi_B * (1-delta_B) * L_B_SS^(1/varphi) * X_B_SS^gamma * gamma;
end
for t = 3*T+1:4*T % x def'n
    A(t,t) = 1;
    A(t,t-3*T) = -C_B_SS/X_B_SS;
    A(t,t-T) = chi_B * delta_B * L_B_SS^(1+1/varphi)/X_B_SS;
end

A_inv = A^(-1);

%----------------------------------------------------------------
% Baseline
%----------------------------------------------------------------

[c_B_base,b_B_base,l_B_base] = c_B_fn(exo,A_inv);

%----------------------------------------------------------------
% Wages
%----------------------------------------------------------------

C_B_w = NaN(T,T);
B_B_w = NaN(T,T);
L_B_w = NaN(T,T);

for t = 1:T
    exo.w_hat(t) = step;
    [c_B_shock,b_B_shock,l_B_shock] = c_B_fn(exo,A_inv);
    C_B_w(:,t) = (c_B_shock-c_B_base)/step;
    B_B_w(:,t) = (b_B_shock-b_B_base)/step;
    L_B_w(:,t) = (l_B_shock-l_B_base)/step;
    exo.w_hat(t) = 0;
end

%----------------------------------------------------------------
% Transfers
%----------------------------------------------------------------

C_B_tau = NaN(T,T);
B_B_tau = NaN(T,T);
L_B_tau = NaN(T,T);

for t = 1:T
    exo.tau_hat(t) = step;
    [c_B_shock,b_B_shock,l_B_shock] = c_B_fn(exo,A_inv);
    C_B_tau(:,t) = (c_B_shock-c_B_base)/step;
    B_B_tau(:,t) = (b_B_shock-b_B_base)/step;
    L_B_tau(:,t) = (l_B_shock-l_B_base)/step;
    exo.tau_hat(t) = 0;
end

%----------------------------------------------------------------
% Nominal Interest Rate
%----------------------------------------------------------------

C_B_ib = NaN(T,T);
B_B_ib = NaN(T,T);
L_B_ib = NaN(T,T);

for t = 1:T
    exo.ib_hat(t) = step;
    [c_B_shock,b_B_shock,l_B_shock] = c_B_fn(exo,A_inv);
    C_B_ib(:,t) = (c_B_shock-c_B_base)/step;
    B_B_ib(:,t) = (b_B_shock-b_B_base)/step;
    L_B_ib(:,t) = (l_B_shock-l_B_base)/step;
    exo.ib_hat(t) = 0;
end

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_B_pi = NaN(T,T);
B_B_pi = NaN(T,T);
L_B_pi = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_B_shock,b_B_shock,l_B_shock] = c_B_fn(exo,A_inv);
    C_B_pi(:,t) = (c_B_shock-c_B_base)/step;
    B_B_pi(:,t) = (b_B_shock-b_B_base)/step;
    L_B_pi(:,t) = (l_B_shock-l_B_base)/step;
    exo.pi_hat(t) = 0;
end

%----------------------------------------------------------------
% Dividends
%----------------------------------------------------------------

C_B_d = NaN(T,T);
B_B_d = NaN(T,T);
L_B_d = NaN(T,T);

for t = 1:T
    exo.d_hat(t) = step;
    [c_B_shock,b_B_shock,l_B_shock] = c_B_fn(exo,A_inv);
    C_B_d(:,t) = (c_B_shock-c_B_base)/step;
    B_B_d(:,t) = (b_B_shock-b_B_base)/step;
    L_B_d(:,t) = (l_B_shock-l_B_base)/step;
    exo.d_hat(t) = 0;
end

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_B_b = NaN(T,T);
B_B_b = NaN(T,T);
L_B_b = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_B_shock,b_B_shock,l_B_shock] = c_B_fn(exo,A_inv);
    C_B_b(:,t) = (c_B_shock-c_B_base)/step;
    B_B_b(:,t) = (b_B_shock-b_B_base)/step;
    L_B_b(:,t) = (l_B_shock-l_B_base)/step;
    exo.zeta_hat(t) = 0;
end

%% DERIVATIVE MATRICES: H TYPES

%----------------------------------------------------------------
% Auxiliary Matrix
%----------------------------------------------------------------

A = zeros(3*T,3*T); % order: (c, l, x)
for t = 1:T % budget constraint
    A(t,t) = C_H_SS;
    A(t,T+t) = -(1-labtax) * W_SS * L_H_SS;
end
for t = T+1:2*T % labor supply relation
    A(t,t) = 1/varphi * chi_H * delta_H * L_H_SS^(1/varphi) + chi_H * (1-delta_H) * L_H_SS^(1/varphi) * X_H_SS^gamma * 1/varphi;
    A(t,t+T) = chi_H * (1-delta_H) * L_H_SS^(1/varphi) * X_H_SS^gamma * gamma;
end
for t = 2*T+1:3*T % x def'n
    A(t,t) = 1;
    A(t,t-2*T) = -C_H_SS/X_H_SS;
    A(t,t-T) = chi_H * delta_H * L_H_SS^(1+1/varphi)/X_H_SS;
end

A_inv = A^(-1);

%----------------------------------------------------------------
% Baseline
%----------------------------------------------------------------

[c_H_base,b_H_base,l_H_base] = c_H_fn(exo,A_inv);

%----------------------------------------------------------------
% Wages
%----------------------------------------------------------------

C_H_w = NaN(T,T);
B_H_w = NaN(T,T);
L_H_w = NaN(T,T);

for t = 1:T
    exo.w_hat(t) = step;
    [c_H_shock,b_H_shock,l_H_shock] = c_H_fn(exo,A_inv);
    C_H_w(:,t) = (c_H_shock-c_H_base)/step;
    B_H_w(:,t) = (b_H_shock-b_H_base)/step;
    L_H_w(:,t) = (l_H_shock-l_H_base)/step;
    exo.w_hat(t) = 0;
end

%----------------------------------------------------------------
% Transfers
%----------------------------------------------------------------

C_H_tau = NaN(T,T);
B_H_tau = NaN(T,T);
L_H_tau = NaN(T,T);

for t = 1:T
    exo.tau_hat(t) = step;
    [c_H_shock,b_H_shock,l_H_shock] = c_H_fn(exo,A_inv);
    C_H_tau(:,t) = (c_H_shock-c_H_base)/step;
    B_H_tau(:,t) = (b_H_shock-b_H_base)/step;
    L_H_tau(:,t) = (l_H_shock-l_H_base)/step;
    exo.tau_hat(t) = 0;
end

%----------------------------------------------------------------
% Nominal Interest Rate
%----------------------------------------------------------------

C_H_ib = NaN(T,T);
B_H_ib = NaN(T,T);
L_H_ib = NaN(T,T);

for t = 1:T
    exo.ib_hat(t) = step;
    [c_H_shock,b_H_shock,l_H_shock] = c_H_fn(exo,A_inv);
    C_H_ib(:,t) = (c_H_shock-c_H_base)/step;
    B_H_ib(:,t) = (b_H_shock-b_H_base)/step;
    L_H_ib(:,t) = (l_H_shock-l_H_base)/step;
    exo.ib_hat(t) = 0;
end

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_H_pi = NaN(T,T);
B_H_pi = NaN(T,T);
L_H_pi = NaN(T,T);

for t = 1:T
    exo.pi_hat(t) = step;
    [c_H_shock,b_H_shock,l_H_shock] = c_H_fn(exo,A_inv);
    C_H_pi(:,t) = (c_H_shock-c_H_base)/step;
    B_H_pi(:,t) = (b_H_shock-b_H_base)/step;
    L_H_pi(:,t) = (l_H_shock-l_H_base)/step;
    exo.pi_hat(t) = 0;
end

%----------------------------------------------------------------
% Dividends
%----------------------------------------------------------------

C_H_d = NaN(T,T);
B_H_d = NaN(T,T);
L_H_d = NaN(T,T);

for t = 1:T
    exo.d_hat(t) = step;
    [c_H_shock,b_H_shock,l_H_shock] = c_H_fn(exo,A_inv);
    C_H_d(:,t) = (c_H_shock-c_H_base)/step;
    B_H_d(:,t) = (b_H_shock-b_H_base)/step;
    L_H_d(:,t) = (l_H_shock-l_H_base)/step;
    exo.d_hat(t) = 0;
end

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_H_b = NaN(T,T);
B_H_b = NaN(T,T);
L_H_b = NaN(T,T);

for t = 1:T
    exo.zeta_hat(t) = step;
    [c_H_shock,b_H_shock,l_H_shock] = c_H_fn(exo,A_inv);
    C_H_b(:,t) = (c_H_shock-c_H_base)/step;
    B_H_b(:,t) = (b_H_shock-b_H_base)/step;
    L_H_b(:,t) = (l_H_shock-l_H_base)/step;
    exo.zeta_hat(t) = 0;
end

%% DERIVATIVE MATRICES: AGGREGATE

%----------------------------------------------------------------
% Wages
%----------------------------------------------------------------

C_w = (1-mu) * C_B_SS/C_SS * C_B_w + mu * C_H_SS/C_SS * C_H_w;
B_w = (1-mu) * B_B_SS/B_SS * B_B_w + mu * B_H_SS/B_SS * B_H_w;
L_w = (1-mu) * L_B_SS/L_SS * L_B_w + mu * L_H_SS/L_SS * L_H_w;

%----------------------------------------------------------------
% Transfers
%----------------------------------------------------------------

C_tau = (1-mu) * C_B_SS/C_SS * C_B_tau + mu * C_H_SS/C_SS * C_H_tau;
B_tau = (1-mu) * B_B_SS/B_SS * B_B_tau + mu * B_H_SS/B_SS * B_H_tau;
L_tau = (1-mu) * L_B_SS/L_SS * L_B_tau + mu * L_H_SS/L_SS * L_H_tau;

%----------------------------------------------------------------
% Nominal Interest Rate
%----------------------------------------------------------------

C_ib = (1-mu) * C_B_SS/C_SS * C_B_ib + mu * C_H_SS/C_SS * C_H_ib;
B_ib = (1-mu) * B_B_SS/B_SS * B_B_ib + mu * B_H_SS/B_SS * B_H_ib;
L_ib = (1-mu) * L_B_SS/L_SS * L_B_ib + mu * L_H_SS/L_SS * L_H_ib;

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

C_pi = (1-mu) * C_B_SS/C_SS * C_B_pi + mu * C_H_SS/C_SS * C_H_pi;
B_pi = (1-mu) * B_B_SS/B_SS * B_B_pi + mu * B_H_SS/B_SS * B_H_pi;
L_pi = (1-mu) * L_B_SS/L_SS * L_B_pi + mu * L_H_SS/L_SS * L_H_pi;

%----------------------------------------------------------------
% Dividends
%----------------------------------------------------------------

C_d = (1-mu) * C_B_SS/C_SS * C_B_d + mu * C_H_SS/C_SS * C_H_d;
B_d = (1-mu) * B_B_SS/B_SS * B_B_d + mu * B_H_SS/B_SS * B_H_d;
L_d = (1-mu) * L_B_SS/L_SS * L_B_d + mu * L_H_SS/L_SS * L_H_d;

%----------------------------------------------------------------
% Demand Shock
%----------------------------------------------------------------

C_b = (1-mu) * C_B_SS/C_SS * C_B_b + mu * C_H_SS/C_SS * C_H_b;
B_b = (1-mu) * B_B_SS/B_SS * B_B_b + mu * B_H_SS/B_SS * B_H_b;
L_b = (1-mu) * L_B_SS/L_SS * L_B_b + mu * L_H_SS/L_SS * L_H_b;

%% DISPLAY RESULTS

C_tau_levels = C_tau * C_SS/Trans_SS;

omega   = C_tau_levels(1,1);
xitheta = C_tau_levels(2,1)/omega;
theta   = (1 + r_b_SS) - xitheta * omega/(1-omega);
xi      = xitheta/theta;

disp(['The average MPC is ' num2str(omega)])
disp(['The average MPE is ' num2str((1-labtax) * W_SS * L_tau(1,1) * L_SS/Trans_SS)])
disp(['The MPC slope is ' num2str(theta)])

disp(['The B MPC is ' num2str(C_B_tau(1,1) * C_B_SS/Trans_B_SS)])
disp(['The B MPE is ' num2str((1-labtax) * W_SS * L_B_tau(1,1) * L_B_SS/Trans_B_SS)])

disp(['The H MPC is ' num2str(C_H_tau(1,1) * C_H_SS/Trans_H_SS)])
disp(['The H MPE is ' num2str((1-labtax) * W_SS * L_H_tau(1,1) * L_H_SS/Trans_H_SS)])

%% SAVE RESULTS

cd([path model task '/_results']);

% aggregate parameters

save param_agg varphi delta_B delta_H chi_B chi_H ...
     epsilon_p phi_p kappa_p ...
     labtax rho_tr phi_pi phi_y phi_dy rho_dr

% household parameters

save param_households beta gamma

% steady state

save SS C_SS L_SS Y_SS Trans_SS G_SS W_SS P_I_SS Pi_SS R_n_SS R_b_SS B_SS D_SS Z_SS r_b_SS

% other quantities

save aux eta alpha mu

% derivative matrices

save jac_matrices C_w B_w L_w C_tau B_tau L_tau C_ib B_ib L_ib C_pi B_pi L_pi C_d B_d L_d C_b B_b L_b

cd([path model task]);