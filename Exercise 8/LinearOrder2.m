function LinearOrder2(X0,A,b)

% A is the interaction matrix size(2,2)
% b is the constant input vector size(2,1)
% X0 is the initial conditions of trajectory

%% input
statesize = 20; % statesize is the size of the illustrated part of statespace
tend = 3; % tend is final time of trajectory evaluation (in multiples of tau)

%%
fs = 14;    % fontsize for legends

%% solve for steady-state
% steady-state condition is 0 = dv/dt = A dot Xss + B or, equivalently, Xss = A^(-1) dot (-B)
Xss = A \ (-b);

%% obtain eigenvectors and eigenvalues
[V, D] = eig( A );

%% alternative variable names
a11=A(1,1);
a12=A(1,2);
a21=A(2,1);
a22=A(2,2);

b1 = b(1);
b2 = b(2);

xss = Xss(1);
yss = Xss(2);

lambda1  = D(1,1)  % eigenvalues
lambda2  = D(2,2)

E1  = V(:,1)       % eigenvectors
E2  = V(:,2)

%% x-range over which nullclines are evaluated (choose symmetrically around initial condition)
xmin = floor(xss - 0.5*statesize);
xmax = ceil(xss + 0.5*statesize);

