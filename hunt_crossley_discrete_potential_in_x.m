function ell = hunt_crossley_discrete_potential_in_x(data, x)
% Computes ell(x; x0) for the H&C model with stiff core.

% Model parameters.
k = data.k;            % Stiffness.
d = data.d;            % Hunt & Crossley dissipation.
xe = data.xe;          % Distance to the "stiff core".
lambda = data.lambda;  % Stiff core exponential parameter.

% Discrete stepping parameters
x0 = data.x0;          % Previous step penetration.
dt = data.dt;          % Time step.

% Potential is constant for x < x0 - xd, i.e. gamma = 0.
% Trick to avoid division by zero:
%   when d = dt/x0, then xd = x0, and xmin = 0.
xd = dt/max(d, dt/x0);  
xmin = max(0, x0 - xd);  % ell = constant for x < x0 - xd.
x = max(xmin, x);

% Dimensionless coordinate.
u = x/xe;

% Dimensionless parameters.
b = 1 - d/dt*x0;
c = d/dt*xe;


factor = k*xe*xe; % Add units.

ell = factor * (((c*u.^2+b*u)*lambda^2+(-2*c*u-b)*lambda+2*c).*exp(u*lambda))/lambda^3;

