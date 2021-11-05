% Course on Spiking Networks
% Teacher: Prof. Dr. Jochen Braun
% Exercise02: Renewal Processes
% Name: Lee Po Shing
% Date: 23/10/2018

clear all
clc
close all

% With size

a = 20; %in ms
S_t_emp = rand(1, 10000); %empirically generated survivor fraction of 1x10000
S_t_emp = sort(S_t_emp);
t_emp = sqrt(-2*a^2*log(S_t_emp)); %intervals generated, in ms
rho_t_emp = t_emp/(a^2); %hazard, in kHz
edges = 0:100; %edges for histogram, in ms
[P_t_emp, edges] = histcounts(t_emp, edges, 'Normalization', 'probability');
bincentre = 0.5:99.5;

S_t_ana = linspace(0, 1, 10001);
S_t_ana = S_t_ana(2:10001); %analytically generated survivor fractions of 1x10000
t_ana = sqrt(-2*a^2*log(S_t_ana)); %intervals generated, in ms
rho_t_ana = t_ana/(a^2); %hazard, in kHz
P_t_ana = rho_t_ana.*exp(-rho_t_ana.*t_ana/2); %interval density, in kHz

figure
axis square
hold on
plot(t_ana, rho_t_ana, 'linewidth', 3, 'linestyle', '--')
plot(t_emp, rho_t_emp)
xlabel('t [ms]')
ylabel('\rho(t) [kHz]')

figure
axis square
hold on
plot(t_ana, P_t_ana, 'linewidth', 3, 'linestyle', '--')
plot(bincentre, P_t_emp)
xlabel('t [ms]')
ylabel('P(t) [kHz]')

figure
axis square
hold on
plot(t_ana, S_t_ana, 'linewidth', 3, 'linestyle', '--')
plot(t_emp, S_t_emp)
xlabel('t [ms]')
ylabel('S(t)')
legend('analytical', 'empirical', 'Location', 'east')

mean_ana = sqrt(pi*a^2/2)
mean_emp = mean(t_emp)
var_ana = 2*a^2-pi*a^2/2
var_emp = var(t_emp)