%% SUFFICIENT STATISTICS APPROXIMATION FOR DIFFERENT VALUES
% Christian Wolf
% this version: 06/19/2024

%% HOUSEKEEPING

clc
clear all
close all

local = '/Users/christianwolf/Dropbox/Research/mpfp_equiv/Codes';
path = [local '/jpe_replication'];
experiment = '/analytical';

addpath(genpath([path '/_auxiliary_functions']));
addpath([path experiment '/_aux']);
addpath([path experiment '/_inputs']);
cd([path experiment]);

%% COMPUTE SUFFICIENT STATISTICS APPROXIMATION

%----------------------------------------------------------------
% Target Demand Paths
%----------------------------------------------------------------

% parameters

T = 500;
r_b_SS = 0.01;

% paths

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

%----------------------------------------------------------------
% Approximation
%----------------------------------------------------------------

% grids

omega_grid = [0.5, 0.4, 0.3, 0.2];
n_omega    = length(omega_grid);

n_theta      = 1;
target_ratio = 1.173758940736347;

tau_equiv_pred = zeros(T,size(demand_target,2),n_omega,n_theta);

% loop

for i_omega = 1:n_omega

    omega      = omega_grid(i_omega);
    theta_grid = target_ratio * (1-omega);
    
    for i_theta = 1:n_theta
    
        theta = theta_grid(i_theta);
        
        xi = (1/omega - 1) * (1 - theta/(1+r_b_SS))/(theta/(1+r_b_SS));
        psi = 1;
        
        C_tau_pred = suffstats_fn(omega,theta,xi,psi,T);
        C_tau_pred_inv = C_tau_pred^(-1);
        
        for i = 1:size(demand_target,2)
            tau_equiv_pred(:,i,i_omega,i_theta) = C_tau_pred^(-1) * demand_target(:,i);
        end
    
    end

end

%% PLOT RESULTS

cd([path experiment '/_results']);

% colors

settings.colors.black  = [0 0 0];
settings.colors.grey   = [230/255 230/255 230/255];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
settings.colors.lnavyblue = 0.25 * [0/255 0/255 50/255] + 0.75 * [1 1 1];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.lpurple = 0.25 * [160/255 32/255 240/255] + 0.75 * [1 1 1];

settings.colors.var = zeros(40,3);
settings.colors.maxhor = 5;
settings.colors.maxweight = 0.8;
for i = 1:40
    if i <= settings.colors.maxhor
        weight = 0 + (i-1) * settings.colors.maxweight/(settings.colors.maxhor-1);
        settings.colors.mpc(i,:) = (1-weight) * settings.colors.orange + weight * [1 1 1];
    else
        settings.colors.mpc(i,:) = (1-weight) * settings.colors.orange + weight * [1 1 1];
    end
end
clear weight

% plot size

plotwidth = 0.25;
gapsize = 0.075;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * (gapsize + plotwidth)];

% figure

max_hor = 21;
scale = 15000/100;

figure(1)

subplot(1,3,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
for i_omega = 1:n_omega
    plot(0:1:max_hor-1,scale*tau_equiv_pred(1:max_hor,1,i_omega),'linewidth',3.5,'linestyle',':','color',settings.colors.mpc(i_omega,:))
    hold on
end
hold on
plot(0:1:max_hor-1,scale*demand_target(1:max_hor,1),'linewidth',5,'linestyle','-','color',settings.colors.grey)
hold on
set(gcf,'color','w')
title('Transitory','interpreter','latex','fontsize',21)
xlabel('Horizon','interpreter','latex','FontSize',17)
ylabel('\$','interpreter','latex','FontSize',17,'Rotation',0)
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
for i_omega = 1:n_omega
    plot(0:1:max_hor-1,scale*tau_equiv_pred(1:max_hor,2,i_omega),'linewidth',3.5,'linestyle',':','color',settings.colors.mpc(i_omega,:))
end
hold on
plot(0:1:max_hor-1,scale*demand_target(1:max_hor,2),'linewidth',5,'linestyle','-','color',settings.colors.grey)
hold on
set(gcf,'color','w')
title('Persistent','interpreter','latex','fontsize',21)
xlabel('Horizon','interpreter','latex','FontSize',17)
% ylabel('\%','interpreter','latex','FontSize',18,'Rotation',0)
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
hold on
plot(0,0,'linewidth',5,'linestyle','-','color',settings.colors.grey)
hold on
for i_omega = 1:n_omega
    plot(0:1:max_hor-1,scale*tau_equiv_pred(1:max_hor,3,i_omega),'linewidth',3.5,'linestyle',':','color',settings.colors.mpc(i_omega,:))
end
hold on
plot(0:1:max_hor-1,scale*demand_target(1:max_hor,3),'linewidth',5,'linestyle','-','color',settings.colors.grey)
hold on
set(gcf,'color','w')
title('Hump-Shaped','interpreter','latex','fontsize',21)
xlabel('Horizon','interpreter','latex','FontSize',17)
% ylabel('\%','interpreter','latex','FontSize',18,'Rotation',0)
legend({'Spending Target','Required Transfer'},'Location','Southeast','fontsize',17,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.25*1.7*pos(3) 1.25*0.7*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('figure_c2','-depsc');

cd([path experiment]);