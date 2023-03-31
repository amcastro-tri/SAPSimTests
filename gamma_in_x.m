function gamma = gamma_in_x(data, x)

% Model parameters.
k = data.k;            % Stiffness.
d = data.d;            % Hunt & Crossley dissipation.
xe = data.xe;          % Distance to the "stiff core".
lambda = data.lambda;  % Stiff core exponential parameter.

% Discrete stepping parameters
x0 = data.x0;          % Previous step penetration.
dt = data.dt;          % Time step.

% Potential is constant for x < x0 - xd, i.e. gamma = 0.
v = (x - x0) / dt;

% Dimensionless coordinate.
u = max(0, x/xe);

factor = k*xe;  % Add units.

bu = u.*exp(lambda*u);

gamma = factor * bu .* max(0, 1 + d * v);
