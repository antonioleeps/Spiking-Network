% Course on Spiking Networks
% Teacher: Prof. Dr. Jochen Braun
% Exercise04: State Space of Deterministic System
% Name: Lee Po Shing
% Date: 25/11/2018

clear all
clc
close all

%% Time-development of random initial states
tau = 10; %time constant
kappa = 0.75; %sigmoidal constant
dt = 1; %time interval
t = 0:dt:150; %time matrix

w = 1; %strength of auto-excitation of E_1

figure
hold on
axis square

for i = 1:10
    [E, E_ss] = timetrajectory(w, t, kappa, dt, tau); %calculation of time-development with a random initial state
    
    plot(E(1,:), E(2,:), 'm', 'linewidth', 1);
    scatter(E(1,1), E(2,1), 'm', 'linewidth', 1.5, 'marker', '*');
    scatter(E_ss(1,:), E_ss(2,:), 'g', 'linewidth', 1.5, 'marker', '+');
end

legend('E', 'E_0', 'E_{ss}', 'location', 'bestoutside')

plotting = 1;
[steadypoint] = plotstatespace(w, kappa, plotting, tau); %calculation and plotting of isoclines, gradient and steady point(s)

%% Effect of w on steady point(s)
w = 1:0.2:4; %array of w
w_temp = [];
E_ss_w = [];

plotting = 0;
for i = 1:length(w)
    [steadypoint] = plotstatespace(w(i), kappa, plotting, tau); %calculation of steady point(s) without plotting
    E_ss_w = [E_ss_w steadypoint];
    w_temp = [w_temp w(i)*ones(1, size(steadypoint, 2))];
end

figure
hold on
scatter(w_temp, E_ss_w(2,:))
scatter(w_temp, E_ss_w(1,:))
legend('E_1', 'E_2')

axis([1 4 -0.2 1.4])
xlabel('w')
ylabel('E')

%% Optional: Jacobian (w = 4)
w = 4;

for j = 1:length(w_temp)
    if w_temp(j) ==w;
        J = zeros(2,2);
        p = w_temp(j)*E_ss_w(2,j) -E_ss_w(1,j) -0.5;
        q = E_ss_w(2,j) -E_ss_w(1,j) -0.5;
        
        J(1,2) = -2*p*kappa^2/(kappa^2 + p^2)^2;
        J(1,1) = -1 - w_temp(j)*J(1,2);
        J(2,1) = 2*q*kappa^2/(kappa^2 + q^2)^2;
        J(2,2) = -1 - J(2,1);
        
        [V, D] = eig(J);
        D
    end
end
    