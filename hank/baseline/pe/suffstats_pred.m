%% ACCURACY OF SUFFICIENT STATISTICS APPROXIMATION
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
model = '/hank/baseline';
task = '/pe';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path model '/_aux']);
addpath([path model '/_inputs/_results']);
cd([path model task]);

%% PREPARATIONS

%----------------------------------------------------------------
% Experiment
%----------------------------------------------------------------

lw_base_indic = 0;
lw_low_indic  = 0;
lw_high_indic = 1;

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

if lw_base_indic == 1

load param_agg
load param_households
load SS
load aux
load jac_matrices

elseif lw_low_indic == 1

load param_agg_low
load param_households_low
load SS_low
load aux_low
load jac_matrices_low

elseif lw_high_indic == 1

load param_agg_high
load param_households_high
load SS_high
load aux_high
load jac_matrices_high

end

%----------------------------------------------------------------
% Time Horizon
%----------------------------------------------------------------

global T

T = size(C_tau,1);

%----------------------------------------------------------------
% Re-definitions
%----------------------------------------------------------------

r_b_SS = R_b_SS - 1;

r_NPV = zeros(T,1);
for t = 1:T
    r_NPV(t) = 1/(1 + r_b_SS)^(t-1);
end

%% CONSTRUCT SUFFICIENT STATISTICS APPROXIMATION

%----------------------------------------------------------------
% Get Sufficient Statistics
%----------------------------------------------------------------

C_tau_true = C_tau * C_SS/Trans_SS;
C_tau_true_inv = C_tau_true^(-1);

omega   = C_tau_true(1,1);
xitheta = C_tau_true(2,1)/omega;
theta   = (1 + r_b_SS) - xitheta * omega/(1-omega);
xi      = xitheta/theta;
psi     = 1;

%----------------------------------------------------------------
% Get Approximation
%----------------------------------------------------------------

C_tau_pred = suffstats_fn(omega,theta,xi,psi,T);
C_tau_pred_inv = C_tau_pred^(-1);

%----------------------------------------------------------------
% Predict Demand Paths
%----------------------------------------------------------------

% get targets

demand_target = NaN(T,3);

demand_target(1,1) = 1;
for t = 2:T
    demand_target(t,1) = 0.3 * demand_target(t-1,1);
end

demand_target(1,2) = 1;
for t = 2:T
    demand_target(t,2) = 0.8 * demand_target(t-1,2);
end

for t = 1:T
    demand_target(t,3) = -exp(-0.1*t) + exp(-0.3*t);
end
demand_target(:,3) = demand_target(:,3) ./ max(abs(demand_target(:,3)));
if demand_target(1,3) < 0
    demand_target(:,3) = -demand_target(:,3);
end

% predict required transfers

tau_equiv_true = zeros(T,size(demand_target,2));
tau_equiv_pred = zeros(T,size(demand_target,2));
for i = 1:size(demand_target,2)
    tau_equiv_true(:,i) = C_tau_true^(-1) * demand_target(:,i);
    tau_equiv_pred(:,i) = C_tau_pred^(-1) * demand_target(:,i);
end

%% PLOT RESULTS

cd([path model task '/_results']);

% colors

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.lnavyblue = 0.25 * [0/255 0/255 50/255] + 0.75 * [1 1 1];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.lpurple = 0.25 * [160/255 32/255 240/255] + 0.75 * [1 1 1];

weights = [1 0.8 0.5 0.2];
settings.colors.orange_all = zeros(length(weights),3);
settings.colors.blue_all = zeros(length(weights),3);
for i = 1:4
    settings.colors.orange_all(i,:) = weights(i) * settings.colors.orange + (1-weights(i)) * [1 1 1];
    settings.colors.blue_all(i,:) = weights(i) * settings.colors.black + (1-weights(i)) * [1 1 1];
end

% plot size

plotwidth = 0.25;
gapsize = 0.075;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * (gapsize + plotwidth)];

% plot entries of derivative matrices

select_hor = [1 6 12 18];

max_hor = 21;

figure(1)
pos = get(gca, 'Position');
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
for i = 1:length(select_hor)
    plot(0:1:max_hor-1,C_tau_true_inv(1:max_hor,select_hor(i)),'linewidth',5,'linestyle','-','color',settings.colors.blue_all(i,:))
    hold on
    plot(0:1:max_hor-1,C_tau_pred_inv(1:max_hor,select_hor(i)),'linewidth',5,'linestyle',':','color',settings.colors.orange_all(i,:))
    hold on
end
set(gcf,'color','w')
xlabel('Horizon','interpreter','latex','FontSize',18)
legend({'Exact HA','Sufficient Statistics'},'Location','Southeast','fontsize',18,'interpreter','latex')
if lw_base_indic == 1
ylim([-5 8])
elseif lw_low_indic == 1
ylim([-1 2])
elseif lw_high_indic == 1
ylim([-15 25])
end
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'auto');
if lw_base_indic == 1
    print('figure_1','-depsc');
elseif lw_low_indic == 1
    print('figure_c3_1','-depsc');
elseif lw_high_indic == 1
    print('figure_c3_2','-depsc');
end

% plot different demand paths

max_hor = 21;
scale = 15000/100;

figure(2)

subplot(1,3,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:max_hor-1,scale*tau_equiv_true(1:max_hor,1),'linewidth',5,'linestyle','-','color',settings.colors.black)
hold on
plot(0:1:max_hor-1,scale*tau_equiv_pred(1:max_hor,1),'linewidth',5,'linestyle',':','color',settings.colors.orange)
hold on
plot(0:1:max_hor-1,scale*demand_target(1:max_hor,1),'linewidth',5,'linestyle','-','color',settings.colors.grey)
hold on
set(gcf,'color','w')
title('Transitory','interpreter','latex','fontsize',23)
xlabel('Horizon','interpreter','latex','FontSize',19)
ylabel('\$','interpreter','latex','FontSize',19,'Rotation',0)
if lw_base_indic == 1
    ylim([-250 750])
    yticks([-250 0 250 500 750])
elseif lw_low_indic == 1
    ylim([-50 250])
    yticks([-50 0 50 100 150 200 250])
elseif lw_high_indic == 1
    ylim([-1000 2500])
    yticks([-1000 -500 0 500 1000 1500 2000 2500])
end
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:max_hor-1,scale*tau_equiv_true(1:max_hor,2),'linewidth',5,'linestyle','-','color',settings.colors.black)
hold on
plot(0:1:max_hor-1,scale*tau_equiv_pred(1:max_hor,2),'linewidth',5,'linestyle',':','color',settings.colors.orange)
hold on
plot(0:1:max_hor-1,scale*demand_target(1:max_hor,2),'linewidth',5,'linestyle','-','color',settings.colors.grey)
hold on
set(gcf,'color','w')
title('Persistent','interpreter','latex','fontsize',23)
xlabel('Horizon','interpreter','latex','FontSize',19)
if lw_base_indic == 1
    ylim([-250 750])
    yticks([-250 0 250 500 750])
elseif lw_low_indic == 1
    ylim([-50 250])
    yticks([-50 0 50 100 150 200 250])
elseif lw_high_indic == 1
    ylim([-1000 2500])
    yticks([-1000 -500 0 500 1000 1500 2000 2500])
end
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0:1:max_hor-1,scale*tau_equiv_true(1:max_hor,3),'linewidth',5,'linestyle','-','color',settings.colors.black)
hold on
plot(0:1:max_hor-1,scale*tau_equiv_pred(1:max_hor,3),'linewidth',5,'linestyle',':','color',settings.colors.orange)
hold on
plot(0:1:max_hor-1,scale*demand_target(1:max_hor,3),'linewidth',5,'linestyle','-','color',settings.colors.grey)
hold on
set(gcf,'color','w')
title('Hump-Shaped','interpreter','latex','fontsize',23)
xlabel('Horizon','interpreter','latex','FontSize',19)
legend({'Exact HA','Sufficient Statistics','Spending Target'},'Location','Northeast','fontsize',19,'interpreter','latex')
if lw_base_indic == 1
    ylim([-250 750])
    yticks([-250 0 250 500 750])
elseif lw_low_indic == 1
    ylim([-50 250])
    yticks([-50 0 50 100 150 200 250])
elseif lw_high_indic == 1
    ylim([-1000 2500])
    yticks([-1000 -500 0 500 1000 1500 2000 2500])
end
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.25*1.7*pos(3) 1.45*0.7*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if lw_base_indic == 1
    print('figure_2','-depsc');
elseif lw_low_indic == 1
    print('figure_c4_1','-depsc');
elseif lw_high_indic == 1
    print('figure_c4_2','-depsc');
end

cd([path model task]);