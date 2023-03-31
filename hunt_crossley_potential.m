function ell = hunt_crossley_potential(data, v)
% Computes ell(v) such that gamma(v) = -dell/dv.
% For small penetrations, the model reduces to Hunt & Crossley.
% For large penetrations, the model penalizes them exponentially.
% See hunt_crossley_force.m
% See symobolic manipulation in hunt_crossley_integrals.wxmx.

k = data.k;        % Stiffness.
d = data.d;            % Hunt & Crossley dissipation.
delta = data.delta;    % Distance to the "stiff core".
lambda = data.lambda;  % Stiff core exponential parameter.
ve = data.ve;          % ve = delta / dt.
vp = data.vp;          % v_phi = -phi0/dt.

vd = 1/d;
v_min = min(vp, vd);

v = min(v_min, v);

ell = (k*delta^2*(((1-d*v)*lambda-d*ve).*((lambda*(vp-v))/ve-1)+d*ve))/(lambda^3).*exp((vp-v)/ve*lambda);
