function t_temp = noisytime(w, t, dt, tau, rho_n, tau_n, E_0, steadypoint)
%Calculation of time-development
A = [w -1;1 -1]; %connection matrix of population E_1 and E_2
b = [-0.5;-0.5]; %external input to both population

%E = zeros(2, length(t)); %matrix of state space at different time values
xi = randn(2,length(t)); %white noise
N_ss = rho_n*sqrt((1-exp(-2*dt/tau_n))/2)*xi; %steady-state of noisy input
N = zeros(2, length(t)); %noisy input to the system

E(:,1) = E_0;

t_temp = 0;

for i = 1:length(t)-1
    N(:,i+1) = N_ss(:,i) + N(:,i)*exp(-dt/tau_n);
    E_ss = Factivation(A*E(:,i) + b) + N(:,i); %calculation of steady-state
    E(:,i+1) = E_ss + (E(:,i) - E_ss)*exp(-dt/tau); %iterative calculation of each state
    
    if ismembertol(E(:,i+1)', steadypoint', 1e-1, 'byrows', true) == [1;1]
        t_temp = t(i+1);
        break
    end
end

% plot(E(1,:), E(2,:))
end