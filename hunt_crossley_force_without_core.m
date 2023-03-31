function gamma = hunt_crossley_force_without_core(data, v)
% gamma(v) according to the H&C model, no stiff core.

k = data.k;            % Stiffness.
d = data.d;            % Hunt & Crossley dissipation.
delta = data.delta;    % Distance to the "stiff core".
ve = data.ve;          % ve = delta / dt.
vp = data.vp;          % v_phi = -phi0/dt.

vd = 1/d;
v_min = min(vp, vd);
v = min(v_min, v);

s = (vp-v)/ve;
bp = s;
impulse_factor = k*delta^2/ve;

gamma = impulse_factor * bp .* (1-d*v); 
