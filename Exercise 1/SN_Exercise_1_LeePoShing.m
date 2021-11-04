% Course on Spiking Networks
% Teacher: Prof. Dr. Jochen Braun
% Exercise01: Random Variables
% Name: Lee Po Shing
% Date: 15/10/2018
% Purpose: Binomial Distribution

clear all
clc
close all

% With size of 20 neurons, the basic properties of the connection matrix,
% including the matrix itself, the two principal eigenvectors and the
% co

N = 3; %number of toss
n = 0:N; %possible total number of heads
p = 0.5; %probability of head

for i = 1:N+1
    f_n(i) = nchoosek(N, n(i))*p^n(i)*(1-p)^(N-n(i)); %binomial distribution
end
F_n = cumsum(f_n); %cummulative binomial distribution

n_r = zeros(1,1000);
for i = 1:N
    n_r = n_r + randi(2, 1, 1000) -1; %empiraical generation of coin-flipping results
end

f_n_e = histcounts(n_r, -0.5:N+0.5, 'Normalization', 'probability'); %
F_n_e = cumsum(f_n_e);

figure
hold on
xlim([0 N+1])
ylim([0 1])
xticks(0:N+1)
plot(n, f_n, 'LineWidth', 2, 'Marker', 'o')
plot(n, F_n, 'LineWidth', 2, 'Marker', '*')

plot(n, f_n_e, 'LineWidth', 1)
plot(n, F_n_e, 'LineWidth', 1)

xlabel('n')
ylabel('P(n)')
legend('f(n)', 'F(n)', 'f(n)_{empirical}', 'F(n)_{empirical}', 'Location', 'East')