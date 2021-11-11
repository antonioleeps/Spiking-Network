function X = Factivation(X)
Fmax = 1;
kappa = 0.75; %sigmoidal constant
%kk = find( X < 0 );
X = Fmax * X.^2 ./ ( kappa^2 + X.^2 );
%X(kk) = zeros(size(kk));

return;