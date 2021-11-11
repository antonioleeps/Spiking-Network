% Course on Spiking Networks
% Teacher: Prof. Dr. Jochen Braun
% Exercise06: Wiener Process
% Name: Lee Po Shing
% Date: 2/1/2019

clear all
clc
close all

tau = 5; %time constant, in ms
x_0 = -70; %in mV
mu = 20; %in mV
rho = 10; %in mV

trial = 1500; %number of realization
step = 80; %number of increments
n_t = randn(trial, step); %samples for white noise
dt = tau/10;
xi_t = n_t / sqrt(dt); %white noise
t = dt*(0:step);

x_t = (mu + rho*sqrt(tau)*xi_t)*dt/tau;
X_t = cumsum(x_t, 2) + x_0;
X_t = [x_0*ones(trial, 1) X_t];

figure
plot(t, X_t)
xlabel('Time (ms)')
ylabel('x(t) (mV)')

%%
t_test = 5;
i = t_test*tau/dt + 1;

figure
hold on
histogram(X_t(:,i), 100, 'Normalization', 'cdf')
xlabel('x (mV)')
ylabel('P(x(5\tau) \leq x)')

w = -30:110; %array of voltage
p_t_w = exp(-(w-mu*t_test-x_0).^2/(2*rho^2*t_test))/sqrt(2*pi*rho^2*t_test); %probability density of voltage at 5*tau
P_t_w = cumsum(p_t_w); %cumulative density

plot(w, P_t_w)
xlim([-40 120])
ylim([0 1.2])
legend('Empirical', 'Theoretical', 'Location', 'Northwest')

%%
omega = 1/(dt/1000) * (0:step/2)/step;

x_t_nodrift = rho*sqrt(tau)*xi_t*dt/tau; %point process
x_t_hat = abs(fft(x_t_nodrift, [], 2)); %FFT of point process

S_omega = x_t_hat.^2/(step*dt); %power spectral density
S_omega = mean(S_omega(1:step/2+1), 1); %neglect negative frequency

S_theo = rho^2/tau;

figure
hold on
plot(omega, S_omega)
plot(omega, mean(S_omega)*ones(size(omega)))
plot(omega, S_theo*ones(size(omega)))
xlabel('Frequency (Hz)')
ylabel('S(\omega)')
legend('Empirical', 'Empirical Avg', 'Theoretical', 'Location', 'Northeast')