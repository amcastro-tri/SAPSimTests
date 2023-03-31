close all
clear all

data.dt = 1.0;               % Time step.
data.x0 = 0.75;                % Initial penetration.

data.k = 1.0;          % Stiffness.
data.d = 4;          % Hunt & Crossley dissipation.
data.xe = 1.0;         % Distance to the "stiff core".
data.lambda = 1;     % Stiff core exponential parameter, dimensionless.

% Print data
data


x = linspace(-0.1, 1.1, 100);

ell = hunt_crossley_discrete_potential_in_x(data, x);
gamma = gamma_in_x(data, x);

x_diff = 0.5*(x(1:end-1)+x(2:end));
gamma_diff = -diff(ell)./diff(x);

figure
plot(x,gamma, x, ell-min(ell), x_diff, gamma_diff,'.');
legend('gamma','ell', '-dell/dx')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse figure with different values of lambda and dissipation.

% lambda > 0, d > 0
data.lambda = 1;
data.d = 4;
vd = 1/data.d;
xd = vd * data.dt;
umin = max(0, data.x0-xd)/data.xe;
u0 = data.x0/data.xe;
gamma_s1 = gamma_in_x(data, x);

% lambda > 0, d = 0
data.lambda = 1;
data.d = 0;
vd = 1/data.d;
xd = vd * data.dt;
gamma_s2 = gamma_in_x(data, x);

% lambda = 0, d > 0
data.lambda = 0;
data.d = 4;
vd = 1/data.d;
xd = vd * data.dt;
gamma_1 = gamma_in_x(data, x);

% lambda = 0, d = 0
data.lambda = 0;
data.d = 0;
vd = 1/data.d;
xd = vd * data.dt;
gamma_2 = gamma_in_x(data, x);


y_min=-1;
y_max=8;
figure
h=plot(x,gamma_s2, 'k-', ...
       x,gamma_s1, 'k--', ...
       x,gamma_2, 'k-', ...
       x,gamma_1, 'k--', ...
    [umin umin],[y_min y_max], 'k-',...
    [u0 u0],[y_min y_max], 'k-');
set(h,'LineWidth',2)
set(h(3:4),'LineWidth',2,'Color',[0.5 0.5 0.5])
set(h(5:6),'LineWidth',1)
axis([min(x) max(x) y_min y_max])

legend(...
'\lambda = 1, x_d/x_e=-\infty', ...
'\lambda = 1, x_d/x_e=0.25',...
'\lambda = 0, x_d/x_e=-\infty', ...
'\lambda = 0, x_d/x_e=0.25',...
'Location','NorthWest');

xlabel('u','FontName','Times', 'FontSize',16)
ylabel('\gamma(u)/(\deltat k x_e)','FontName','Times', 'FontSize',16)
set(gca,'FontName','Times','FontSize',16)

% Second set of axes for additional ticks
ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
set(ax2,'FontName','Times','FontSize',16)

xticks(ax2,[0.5 0.75])
xticklabels(ax2,{'u_0-u_d', 'u_0'})
yticks(ax2,[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potential figure with different values of lambda and dissipation.
% lambda > 0, d > 0
data.lambda = 1;
data.d = 4;
vd = 1/data.d;
xd = vd * data.dt;
umin = max(0, data.x0-xd)/data.xe;
u0 = data.x0/data.xe;
gamma_s1 = hunt_crossley_discrete_potential_in_x(data, x);
gamma_s1 = gamma_s1-gamma_s1(1);

% lambda > 0, d = 0
data.lambda = 1;
data.d = 0;
vd = 1/data.d;
xd = vd * data.dt;
gamma_s2 = hunt_crossley_discrete_potential_in_x(data, x);
gamma_s2 = gamma_s2-gamma_s2(1);

% lambda = 0, d > 0
data.lambda = 1e-3;
data.d = 4;
vd = 1/data.d;
xd = vd * data.dt;
gamma_1 = hunt_crossley_discrete_potential_in_x(data, x);
gamma_1 = gamma_1-gamma_1(1);

% lambda = 0, d = 0
data.lambda = 1e-3;
data.d = 0;
vd = 1/data.d;
xd = vd * data.dt;
gamma_2 = hunt_crossley_discrete_potential_in_x(data, x);
gamma_2 = gamma_2-gamma_2(1);

y_min=-0.25;
y_max=2;
figure
h=plot(x,gamma_s2, 'k-', ...
       x,gamma_s1, 'k--', ...
       x,gamma_2, 'k-', ...
       x,gamma_1, 'k--', ...
    [umin umin],[y_min y_max], 'k-',...
    [u0 u0],[y_min y_max], 'k-');
set(h,'LineWidth',2)
set(h(3:4),'LineWidth',2,'Color',[0.5 0.5 0.5])
set(h(5:6),'LineWidth',1)
axis([min(x) max(x) y_min y_max])

legend(...
'\lambda = 1, x_d/x_e=-\infty', ...
'\lambda = 1, x_d/x_e=0.25',...
'\lambda = 0, x_d/x_e=-\infty', ...
'\lambda = 0, x_d/x_e=0.25',...
'Location','NorthWest');

xlabel('u','FontName','Times', 'FontSize',16)
ylabel('B(u) = ℓ(u)/(k x_e^2)','FontName','Times', 'FontSize',16)
set(gca,'FontName','Times','FontSize',16)

% Second set of axes for additional ticks
ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
% set the same Limits and Ticks on ax2 as on ax1;
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'));
set(ax2,'FontName','Times','FontSize',16)

xticks(ax2,[0.5 0.75])
xticklabels(ax2,{'u_0-u_d', 'u_0'})
yticks(ax2,[])


