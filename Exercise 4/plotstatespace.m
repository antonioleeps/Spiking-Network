function [steadypoint] = plotstatespace(w, kappa, plotting, tau) 
%% Isocline for dE_1/dt = 0
E1 = 0:0.02:0.98;
E_2_plus = w*E1-0.5 + 3*sqrt(-E1./(E1-1))/4;
E_2_minus = w*E1-0.5 - 3*sqrt(-E1./(E1-1))/4;

%% Isocline for dE_2/dt = 0
E2 = 0:0.02:0.98;
E_1_plus = E2+0.5 + 3*sqrt(-E2./(E2-1))/4;
E_1_minus = E2+0.5 - 3*sqrt(-E2./(E2-1))/4;

%% Steady-point
steadypoint = [];

dE_1_E2 = -E_1_plus + Factivation((w*E_1_plus - E2 - 0.5), 1, kappa);
idk = [];
for i = 1:length(dE_1_E2)-1
    if dE_1_E2(i)*dE_1_E2(i+1) <0
        idk = [idk i];
    end
end

for i = 1:length(idk)
    E2_temp = E2(idk(i)):0.001:E2(idk(i)+1);
    E_1_temp = E2_temp+0.5 + 3*sqrt(-E2_temp./(E2_temp-1))/4;
    dE_1_E2_temp = -E_1_temp + Factivation((w*E_1_temp - E2_temp - 0.5), 1, kappa);
    
    j = 1;
    while dE_1_E2_temp(j)*dE_1_E2_temp(j+1) >=0
        j = j+1;
    end
    
    E2_temp_2 = (E2_temp(j) + E2_temp(j+1))/2;
    E_1_plus_2 = E2_temp_2 + 0.5 + 3*sqrt(-E2_temp_2./(E2_temp_2-1))/4;
    
    steadypoint = [steadypoint [E2_temp_2;E_1_plus_2]];
end

dE_1_E2 = -E_1_minus + Factivation((w*E_1_minus - E2 - 0.5), 1, kappa);
idk = [];
for i = 1:length(dE_1_E2)-1
    if dE_1_E2(i)*dE_1_E2(i+1) <0
        idk = [idk i];
    end
end

for i = 1:length(idk)
    E2_temp = E2(idk(i)):0.001:E2(idk(i)+1);
    E_1_temp = E2_temp+0.5 - 3*sqrt(-E2_temp./(E2_temp-1))/4;
    dE_1_E2_temp = -E_1_temp + Factivation((w*E_1_temp - E2_temp - 0.5), 1, kappa);
    
    j = 1;
    while dE_1_E2_temp(j)*dE_1_E2_temp(j+1) >=0
        j = j+1;
    end
    
    E2_temp_2 = (E2_temp(j) + E2_temp(j+1))/2;
    E_1_plus_2 = E2_temp_2 + 0.5 - 3*sqrt(-E2_temp_2./(E2_temp_2-1))/4;
    
    steadypoint = [steadypoint [E2_temp_2;E_1_plus_2]];
end

%% Gradient
while plotting == 1
    xv = linspace(-0.2, 1.4, 33);
    yv = flip(xv);
    
    [XV, YV] = meshgrid(xv, yv);
    
    F1 = (-XV + Factivation((w*XV - YV - 0.5), 1, kappa))/tau;
    F2 = (-YV + Factivation((XV - YV - 0.5), 1, kappa))/tau;
    
    %% Plotting of the properties of state-space
    quiver(XV, YV, F1, F2, 'k', 'DisplayName', 'Gradient'); %gradients
    plot(E1, E_2_plus, 'r', 'DisplayName', 'dE_1/dt = 0');
    plot(E1, E_2_minus, 'r', 'HandleVisibility', 'off');
    plot(E_1_plus, E2, 'b', 'DisplayName', 'dE_2/dt = 0');
    plot(E_1_minus, E2, 'b', 'HandleVisibility','off');
    
    scatter(steadypoint(2,:), steadypoint(1,:), 'o', 'DisplayName', 'E_{ss}');
    
    axis([-0.2 1.4 -0.2 1.4])
    xlabel('E_1')
    ylabel('E_2')
    
    plotting = 0;
end
end