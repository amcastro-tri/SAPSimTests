close all

dt = 0.01;               % Time step.
x0 = 5e-4;               % Initial penetration.

data.k = 100.0;          % Stiffness.
data.d = 0.1;            % Hunt & Crossley dissipation.
data.delta = 1.0e-3;      % Distance to the "stiff core".
data.lambda = 2.0;       % Stiff core exponential parameter, dimensionless.
data.ve = data.delta/dt; % ve = delta / dt.
data.vp = x0 / dt;       % v_phi = -phi0/dt.
data.vd = 1/data.d;

% Print data
data

vd = 1/data.d;
v_min = min(data.vp, vd);

%v = linspace(-1.5, 0, 1000)*data.ve + v_min;
v = linspace(-1.5, 0, 1000)*data.vp + v_min;

ell = hunt_crossley_potential(data, v);
gamma = hunt_crossley_force(data, v);
gamma_no_core = hunt_crossley_force_without_core(data, v);

figure
plot(v,gamma, v, gamma_no_core);

% Dimensionless
impulse_scale = data.k*data.delta^2/data.ve;
energy_scale = data.k*data.delta^2;
gamma = gamma / impulse_scale;
gamma_no_core = gamma_no_core / impulse_scale;
ell = ell / energy_scale;
v = (v-v_min)/data.ve;

dell_dv = diff(ell)/(v(2)-v(1));

figure
plot(v,gamma, v(2:end), -dell_dv, v, gamma_no_core);
legend('w/core','-dl/dv', 'no/core')


