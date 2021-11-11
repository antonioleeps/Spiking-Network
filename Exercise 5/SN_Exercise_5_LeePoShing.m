% Course on Spiking Networks
% Teacher: Prof. Dr. Jochen Braun
% Exercise05: Autocorrelation and Power
% Name: Lee Po Shing
% Date: 21/12/2018

clear all
clc
close all

nu = 100; %spiking rate, in Hz
a = sqrt(2/pi) / nu; %in s

S_t = rand(1, 200); %empirically generated survivor fraction of 1x10000
t_isi = sqrt(-2*a^2*log(S_t)); %intervals generated, in s
t_i = round(cumsum(t_isi), 4); %spike time, in s

dt = 1e-4; %discrete interval, in s
alpha = nu*dt; %spike in dt
t_end = t_i(end); %length of signal, in s

t = round(0:dt:t_end, 4); %time array for realization
x_i = ismember(t, t_i) - alpha; %step interval
X_t = cumsum(x_i); %realization

x_i = x_i - mean(x_i); %zero-mean
i = length(t);
x_hat = fft(x_i, i);
omega = 1/dt *(0:i/2)/i ;

S_omega = x_hat.*conj(x_hat) / t_end;
S_omega = S_omega(1:i/2 +1);

figure
plot(log10(omega), S_omega)
xlabel('log_{10} Frequency (Hz)')
ylabel('S(\omega)')

figure
plot(t, X_t);
xlabel('Time (s)')
ylabel('X(t)')

%%
tau = 0:dt:0.04; %time vector, in s
f_tau = (tau/ a^2).*exp(-tau.^2 /(2*a^2)); %interval distribution
%f_tau(tau < 0) = 0;

j = length(tau);
f_hat = fft(f_tau, j)/sqrt(j);

C_hat = nu * real((1+f_hat)./(1-f_hat));
C_hat = C_hat(1:j/2+1);
C_omega = 1/dt * (0:j/2)/j;

figure
plot(log10(C_omega), C_hat)
xlabel('log_{10} Frequency (Hz)')
ylabel('Power')

figure
plot(tau, f_tau);
xlabel('\tau (s)')
ylabel('f(\tau)')

%v = nu^2;
%v_hat = fft(v, j);
%R_hat = C_hat + 2*v_hat(1:j/2+1);

%figure
%plot(log10(C_omega), R_hat)