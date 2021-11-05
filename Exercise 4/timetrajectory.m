function [E, E_ss] = timetrajectory(w, t, kappa, dt, tau)
%Calculation of time-development of a randomly generated initial state

A = [w -1;1 -1]; %connection matrix of population E_1 and E_2
b = [-0.5;-0.5]; %external input to both population

E = zeros(2, length(t)); %vector of state space at different time values
E_0 = 1.6*rand(2, 1)-0.2; %random generation of the initial state

E(:,1) = E_0;

for i = 1:length(t)-1
    E_ss = Factivation((A*E(:,i) + b), 1, kappa); %calculation of steady-state
    E(:,i+1) = E_ss + (E(:,i) - E_ss)*exp(-dt/tau); %iterative calculation of each state
end

end