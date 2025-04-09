radius_core = 1.0;
radius_total = 1.1;

rho_core_fluid = 998.2;
rho_core_solid = 2100;
cp_core_fluid = 4182;
cp_core_solid = 837;
 
porosity = 0.05;
rho_core = porosity * rho_core_fluid + (1 - porosity) * rho_core_solid;
V_core = 4/3*pi*radius_core^3;
m_core = rho_core*V_core;
cp_core = (porosity * cp_core_fluid * rho_core_fluid + (1 - porosity) * cp_core_solid * rho_core_solid) / rho_core;

rho_insulation = 44;
V_ins = 4/3*pi*radius_total^3-V_core;
m_ins = rho_insulation*V_ins;
cp_ins = 700;

% Define constants
oneday = 24*60*60;
cp_f = 4182;        % Fluid specific heat capacity (J/kg·K)
Tamb = 263.15;        % Ambient temperature (K)
Tfi = @(t) (t <= 60*oneday) * (90 + 273.15) + (t > 60*oneday & t <= 80*oneday) * (11+ 273.15) + (t > 80*oneday) * (5 + 273.15);
mdot_f = @(t) (t <= 60*oneday) * 4 + (t > 60*oneday & t <= 80*oneday) * 0.0 + (t > 80*oneday) * 6;

% Thermal resistances
R_hexc = 0.231;
R_ci = 0.152;
R_amb = 0.107;

% Define the DAE system
function dydt = dae_system(t, y, params)
    % Extract constants from params
    m_core = params.m_core;
    cp_core = params.cp_core;
    m_ins = params.m_ins;
    cp_ins = params.cp_ins;
    R_hexc = params.R_hexc;
    R_ci = params.R_ci;
    R_amb = params.R_amb;
    mdot_f = params.mdot_f;
    cp_f = params.cp_f;
    Tamb = params.Tamb;
    Tfi = params.Tfi;

    % Unpack variables
    Tcore = y(1);
    Tins = y(2);
    Tfo = y(3);
    
    % ODEs
    dTcore_dt = -(Tcore - Tfo) / R_hexc - (Tcore - Tins) / R_ci;
    dTins_dt =  -(Tins - Tamb) / R_amb + (Tcore - Tins) / R_ci;
    
    % Return the derivatives
    dydt = zeros(3,1);
    dydt(1) = dTcore_dt / (m_core * cp_core);
    dydt(2) = dTins_dt / (m_ins * cp_ins);
    dydt(3) = (Tfo-Tcore)/R_hexc - mdot_f(t)*cp_f*(Tfi(t)-Tfo); 
end

% Define parameters struct to pass to the DAE system
params = struct('m_core', m_core, 'cp_core', cp_core, ...
                'm_ins', m_ins, 'cp_ins', cp_ins, ...
                'R_hexc', R_hexc, 'R_ci', R_ci, 'R_amb', R_amb, ...
                'mdot_f', mdot_f, 'cp_f', cp_f, ...
                'Tamb', Tamb, 'Tfi', Tfi);

%options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
M = [1 0 0; 0 1 0; 0 0 0];
options = odeset('Mass',M,'RelTol', 1e-6, 'AbsTol', 1e-8);

% ------------------------------------ CHARGE ---------------------------

% Initial conditions
Tcore0 = 293.15;      % Initial temperature of the core (K)
Tins0 = 293.15;       % Initial temperature of the insulation (K)
Tfo0 = (Tcore0+R_hexc*Tfi(0)*cp_f*mdot_f(0))/(R_hexc*cp_f*mdot_f(0) + 1);

% Define time span
tspan = [0, 60*oneday]; % Time span for the simulation (seconds)
y0 = [Tcore0; Tins0; Tfo0];
[t, y] = ode15s(@(t, y) dae_system(t, y, params), tspan, y0, options);
t_all = t;

% Extract solutions
Tcore = y(:,1);
Tins = y(:,2);
Tfo = y(:,3);


% ------------------------------------- STORAGE -----------------------

% Initial conditions
tstart = 60*oneday+0.001;
Tcore0 = Tcore(end);      % Initial temperature of the core (K)
Tins0 = Tins(end);       % Initial temperature of the insulation (K)
Tfo0 = (Tcore0+R_hexc*Tfi(tstart)*cp_f*mdot_f(tstart))/(R_hexc*cp_f*mdot_f(tstart) + 1);

tspan = [tstart, 60*oneday+20*oneday]; % Time span for the simulation (seconds)
y0 = [Tcore0; Tins0; Tfo0];
[t, y] = ode15s(@(t, y) dae_system(t, y, params), tspan, y0, options);
t_all = [t_all; t];

% Extract solutions
Tcore = [Tcore; y(:,1)];
Tins = [Tins; y(:,2)];
Tfo = [Tfo; y(:,1)]; % Tfo = Tcore during storage


% ----------------------------------- DISCHARGE -----------------------

% Initial conditions
tstart = 80*oneday+0.001;
Tcore0 = Tcore(end);      % Initial temperature of the core (K)
Tins0 = Tins(end);       % Initial temperature of the insulation (K)
Tfo0 = (Tcore0+R_hexc*Tfi(tstart)*cp_f*mdot_f(tstart))/(R_hexc*cp_f*mdot_f(tstart) + 1);

tspan = [tstart, 110*oneday]; % Time span for the simulation (seconds)
y0 = [Tcore0; Tins0; Tfo0];
[t, y] = ode15s(@(t, y) dae_system(t, y, params), tspan, y0, options);
t_all = [t_all; t];

% Extract solutions
Tcore = [Tcore; y(:,1)];
Tins = [Tins; y(:,2)];
Tfo = [Tfo; y(:,3)];

% ------------------------------------- VISUALIZE ---------------------
T_2D = importfile('zone-temperatures-rfile-axisymmetric.out');

% Plot the results
figure;

% Set figure size (15 cm x 15 cm)
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 20, 20]);

plot(t_all/oneday, Tcore-273.15, T_2D.time/oneday, T_2D.Core-273.15, t_all/oneday, Tins-273.15, T_2D.time/oneday, T_2D.Insulation-273.15,'LineWidth',2);
grid on
set(gca,'ylim',[-10,45])
set(gca,'xlim',[0,95])

% Add labels and legend
xlabel('Time [days]', 'FontSize', 14); % Increase font size for x-axis label
ylabel('Temperature [^\circC]', 'FontSize', 14); % Increase font size for y-axis label
legend({'Core average (DAE)', 'Core (CFD)', 'Insulation average (DAE)', 'Insulation (CFD)'}, ...
       'FontSize', 12, 'Location', 'southwest'); % Set legend to southwest


