% Course on Spiking Networks
% Teacher: Prof. Dr. Jochen Braun
% Exercise08: State Space of Stochastic System
% Name: Lee Po Shing
% Date: 17/01/2019

clear all
clc
close all

%% Time-development of random initial states
tau = 10; %time constant
tau_n = 2;
w = 3.6; %strength of auto-excitation of E_1
rho_n = 1;
dt = 0.05*tau_n; %time interval

% t_end = 20*tau;
% t = 0:dt:t_end; %time matrix
% E_0 = [0.1;0.1];
% trial = 20;
% 
% figure
% hold on
% axis square
% 
% for i = 1:trial
%     [N, E] = noisytimetrajectory(w, t, dt, tau, rho_n, tau_n, E_0); %calculation of time-development with a random initial state
%     plot(E(1,:), E(2,:), 'linewidth', 0.1);
%     scatter(E(1,1), E(2,1), 'm', 'linewidth', 0.5, 'marker', '*');
%     scatter(E(1,end), E(2,end), 'g', 'linewidth', 0.5, 'marker', '+');
% end
% 
% legend('E', 'E_0', 'E_{end}', 'location', 'bestoutside')
% plotting = 1;
% [steadypoint] = plotstatespace(w, plotting, tau); %calculation and plotting of isoclines, gradient and steady point(s)
% 
% %% Shorter time, Higher noise
% rho_n = 2;
% t_end = 10*tau;
% t = 0:dt:t_end; %time matrix
% trial = 150;
% 
% E_0 = [0.9;0.1];
% t_temp = zeros(1, trial);
% 
% for i = 1:trial
%      t_temp(i) = noisytime(w, t, dt, tau, rho_n, tau_n, E_0, steadypoint(:,3)); %calculation of time-development with a random initial state
% end
% 
% figure
% hold on
% histogram(t_temp(find(t_temp)), 100, 'Normalization', 'cdf');
% 
% E_0 = [0.2;0.6];
% t_temp = zeros(1, trial);
% 
% for i = 1:trial
%      t_temp(i) = noisytime(w, t, dt, tau, rho_n, tau_n, E_0, steadypoint(:,1)); %calculation of time-development with a random initial state
% end
% 
% histogram(t_temp(find(t_temp)), 100, 'Normalization', 'cdf');
% xlabel('Time')
% ylabel('Cumulative distribution')
% legend('High to low', 'Low to high')
% 
% %% Longer time, lower noise
% rho_n = 1;
% t_end = 50*tau;
% t = 0:dt:t_end; %time matrix
% 
% E_0 = [0.9;0.1];
% t_temp = zeros(1, trial);
% 
% for i = 1:trial
%      t_temp(i) = noisytime(w, t, dt, tau, rho_n, tau_n, E_0, steadypoint(:,3)); %calculation of time-development with a random initial state
% end
% 
% figure
% hold on
% histogram(t_temp(find(t_temp)), 100, 'Normalization', 'cdf');
% 
% E_0 = [0.2;0.6];
% t_temp = zeros(1, trial);
% 
% for i = 1:trial
%      t_temp(i) = noisytime(w, t, dt, tau, rho_n, tau_n, E_0, steadypoint(:,1)); %calculation of time-development with a random initial state
% end
% 
% histogram(t_temp(find(t_temp)), 100, 'Normalization', 'cdf');
% xlabel('Time')
% ylabel('Cumulative distribution')
% legend('High to low', 'Low to high')

%% Spectral analysis
step = 100;
t_end = step*dt;
t = 0:dt:t_end; %time vector
rho_n = 0.1;
trial = 125;

E_0 = [0.9;0.1];
E_all = zeros(2, length(t), trial);
N_all = E_all;

for i = 1:trial
    [N_all(:,:,i), E_all(:,:,i)] = noisytimetrajectory(w, t, dt, tau, rho_n, tau_n, E_0); %calculation of time-development with a random initial state
end

omega = 1/dt * (0:step/2)/step;

E_hat = abs(fft(E_all, [], 2));
S_omega = mean(E_hat.^2/t_end, 3); %power spectral density
S_omega_E_high = S_omega(:, 1:step/2+1); %neglect negative frequency

N_hat = abs(fft(N_all, [], 2));
S_omega = mean(N_hat.^2/t_end, 3); %power spectral density
S_omega_N = S_omega(:, 1:step/2+1); %neglect negative frequency

figure
hold on
plot(omega, S_omega_N)
y = tau_n*rho_n^2./(1+tau_n^2*omega.^2);
plot(omega, y)
xlabel('\omega')
ylabel('S(\omega)')
legend('N_1', 'N_2', 'Theoretical')

E_0 = [0.2;0.6];
E_all = zeros(2, length(t), trial);

for i = 1:trial
    [N_all(:,:,i), E_all(:,:,i)] = noisytimetrajectory(w, t, dt, tau, rho_n, tau_n, E_0); %calculation of time-development with a random initial state
end

E_hat = abs(fft(E_all, [], 2));
S_omega = mean(E_hat.^2/t_end, 3); %power spectral density
S_omega_E_low = S_omega(:, 1:step/2+1); %neglect negative frequency

figure
hold on
plot(log(omega), S_omega_E_low-S_omega_E_high)
xlabel('log \omega')
ylabel('\Delta S(\omega)_{low-high}')