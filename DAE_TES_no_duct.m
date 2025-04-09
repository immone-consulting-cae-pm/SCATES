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
cp_f = 4182;        % Fluid specific heat capacity (J/kg·K)
Tamb = 263.15;        % Ambient temperature (K)
Tfi = @(t) (t <= 180*60) * (90 + 273.15) + (t > 60*180 & t <= 60*900) * (20 + 273.15) + (t > 60*900) * (5 + 273.15);
mdot_f = @(t) (t <= 180*60) * 1 + (t > 60*180 & t <= 60*900) * 0 + (t > 60*900) * 3;

% Thermal resistances
R_hexc = 1e-7;
R_ci = 0.09;
R_amb = 0.08;

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
    Tcontrol = 0.5*(Tfi(t)+Tfo);
    
    % ODEs
    dTcore_dt = -(Tcore - Tcontrol) / R_hexc - (Tcore - Tins) / R_ci;
    dTins_dt =  -(Tins - Tamb) / R_amb + (Tcore - Tins) / R_ci;
    
    % Return the derivatives
    dydt = zeros(3,1);
    dydt(1) = dTcore_dt / (m_core * cp_core);
    dydt(2) = dTins_dt / (m_ins * cp_ins);
    dydt(3) = (Tcontrol-Tcore)/R_hexc - mdot_f(t)*cp_f*(Tfi(t)-Tfo); 
end

% Define parameters struct to pass to the DAE system
params = struct('m_core', m_core, 'cp_core', cp_core, ...
                'm_ins', m_ins, 'cp_ins', cp_ins, ...
                'R_hexc', R_hexc, 'R_ci', R_ci, 'R_amb', R_amb, ...
                'mdot_f', mdot_f, 'cp_f', cp_f, ...
                'Tamb', Tamb, 'Tfi', Tfi);

M = [1 0 0; 0 1 0; 0 0 0];
options = odeset('Mass',M,'RelTol', 1e-6, 'AbsTol', 1e-8);

% ------------------------------------ CHARGE ---------------------------

% Initial conditions
Tcore0 = 293.15;      % Initial temperature of the core (K)
Tins0 = 293.15;       % Initial temperature of the insulation (K)
Tfo0 = (Tcore0+R_hexc*Tfi(0)*cp_f*mdot_f(0))/(R_hexc*cp_f*mdot_f(0) + 1);

y0 = [Tcore0; Tins0; Tfo0];

% Define time span
tspan = [0, 3*3600]; % Time span for the simulation (seconds)

[t, y] = ode15s(@(t, y) dae_system(t, y, params), tspan, y0, options);
t_all = t;

% Extract solutions
Tcore = y(:,1);
Tins = y(:,2);
Tfo = y(:,3);


% ------------------------------------- STORAGE -----------------------

% Initial conditions
Tcore0 = Tcore(end);      % Initial temperature of the core (K)
Tins0 = Tins(end);       % Initial temperature of the insulation (K)
Tfo0 = Tcore0;
tspan = [3*3600+0.001, 15*3600]; % Time span for the simulation (seconds)

y0 = [Tcore0; Tins0; Tfo0];

[t, y] = ode15s(@(t, y) dae_system(t, y, params), tspan, y0, options);
t_all = [t_all; t];

% Extract solutions
Tcore = [Tcore; y(:,1)];
Tins = [Tins; y(:,2)];
Tfo = [Tfo; y(:,3)];


% ----------------------------------- DISCHARGE -----------------------

% Initial conditions
tstart = 15*3600+0.001;
Tcore0 = Tcore(end);      % Initial temperature of the core (K)
Tins0 = Tins(end);       % Initial temperature of the insulation (K)
Tfo0 = (Tcore0+R_hexc*Tfi(tstart)*cp_f*mdot_f(tstart))/(R_hexc*cp_f*mdot_f(tstart) + 1);
tspan = [tstart, 17*3600]; % Time span for the simulation (seconds)

y0 = [Tcore0; Tins0; Tfo0];

[t, y] = ode15s(@(t, y) dae_system(t, y, params), tspan, y0, options);
t_all = [t_all; t];

% Extract solutions
Tcore = [Tcore; y(:,1)];
Tins = [Tins; y(:,2)];
Tfo = [Tfo; y(:,3)];

% ------------------------------------- VISUALIZE ---------------------
%T = importfile('laminarzone-zone-temperatures-rfile.out');
T = importfile('turbulent_zone-temperatures-rfile.out');

% Plot the results
figure;

% Set figure size (15 cm x 15 cm)
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 20, 20]);

plot(t_all/3600, Tcore-273.15, T.time/3600, T.Core-273.15, t_all/3600, Tins-273.15, T.time/3600, T.Insulation-273.15,'LineWidth',2);
grid on
set(gca,'xlim',[0,17])

% Add labels and legend
xlabel('Time [h]', 'FontSize', 14); % Increase font size for x-axis label
ylabel('Temperature [^\circC]', 'FontSize', 14); % Increase font size for y-axis label
legend({'Core average (DAE)', 'Core (CFD)', 'Insulation average (DAE)', 'Insulation (CFD)'}, ...
       'FontSize', 12, 'Location', 'southwest'); % Set legend to southwest


