clc; clear; close all;

%% Balloon and Payload Specifications
ascent_speed = 5;            % Balloon ascent speed (m/s)
descent_speed = 4;           % Balloon descent speed (m/s)
payload_mass = 3;            % Payload mass (kg)
payload_specific_heat = 1300;% Specific heat capacity (J/kgK)
payload_volume = 0.003;      % Payload volume (m³)

%% Thermal Properties
thermal_constants = struct(...
    'stefan_boltzmann', 5.670e-8, ...
    'air_thermal_conductivity', 0.025, ...
    'foam_emissivity', 0.6, ...
    'foam_absorptivity', 0.15, ...
    'mylar_emissivity', 0.05, ...
    'mylar_absorptivity', 0.1 ...
);

%% Geometry
payload_geometry = struct(...
    'length', 0.3, ...
    'width', 0.1, ...
    'height', 0.1 ...
);

% Calculate payload surface area
payload_surface_area = 2 * (payload_geometry.length * payload_geometry.width + ...
                             payload_geometry.length * payload_geometry.height + ...
                             payload_geometry.width * payload_geometry.height);

%% Environmental Conditions
solar_constants = struct(...
    'solar_flux', 1361, ...
    'albedo', 0.3, ...
    'visibility_factor', 0.15, ...
    'planetary_flux', 220 ...
);

%% Simulation Parameters
max_altitude = 38000;        % Maximum altitude (m)
simulation_time = 50 + (max_altitude / ascent_speed + max_altitude / descent_speed);  % Total round trip time
time_step = 0.1;             % Time step for calculations (seconds)
target_temp_min = 0;         % Minimum target temperature (°C)
target_temp_max = 10;        % Maximum target temperature (°C)

%% Initialize Arrays
num_steps = floor(simulation_time / time_step);
time = linspace(0, simulation_time, num_steps);

% Create altitude profile with ascent and descent
altitude = zeros(size(time));
for i = 1:length(time)
    % Ascent phase
    if time(i) <= max_altitude / ascent_speed
        altitude(i) = ascent_speed * time(i);
    % Descent phase
    elseif time(i) <= max_altitude / ascent_speed + max_altitude / descent_speed
        altitude(i) = max_altitude - descent_speed * (time(i) - max_altitude / ascent_speed);
    else
        altitude(i) = 0;
    end
end

% Initialize temperature arrays
temperature_external = zeros(size(time));
temperature_payload = zeros(size(time));
heater_power = zeros(size(time));

%% Calculate External Temperature Profile
for i = 1:length(time)
    current_altitude = altitude(i);

    if current_altitude < 11000
        temperature_external(i) = 288.15 - (0.0065 * current_altitude); % No conversion to Celsius yet
    elseif current_altitude < 20000
        temperature_external(i) = 216.65; % Constant temperature
    else
        temperature_external(i) = 216.65 + (current_altitude - 20000) * 0.001; % Linear increase
    end
end

temperature_external = temperature_external - 273.15; % Convert to Celsius after calculation

%% Initial Conditions
temperature_payload(1) = temperature_external(1);

%% Thermal Simulation
for i = 2:length(time)
    % Heat Transfer Calculations
    Q_solar = thermal_constants.foam_absorptivity * solar_constants.solar_flux * payload_surface_area / 2;
    Q_albedo = thermal_constants.foam_absorptivity * solar_constants.solar_flux * solar_constants.albedo * payload_surface_area / 2;
    Q_planetary = thermal_constants.foam_emissivity * solar_constants.planetary_flux * payload_surface_area / 2;
    
    % Radiation Heat Loss
    Q_radiation = thermal_constants.foam_emissivity * thermal_constants.stefan_boltzmann * ...
                  payload_surface_area * (temperature_payload(i-1) + 273.15)^4;
    
    % Convection Heat Loss (simplified model)
    if altitude(i) < 10000
        h_convection = thermal_constants.air_thermal_conductivity / payload_geometry.length;
        Q_convection = h_convection * (temperature_payload(i-1) - temperature_external(i-1));
    else
        Q_convection = 0;
    end
    
    %Heater Power Control
    if temperature_payload(i-1) < target_temp_min
        heater_power(i) = 5;  % Maximum power if below minimum
    elseif temperature_payload(i-1) > target_temp_max
        heater_power(i) = 0;   % No heating if above maximum
    else
        %Proportional control
        heater_power(i) = 5 * (target_temp_max - temperature_payload(i-1)) / target_temp_max;
    end

    
    % Net Heat Balance
    Q_net = Q_solar + Q_albedo + Q_planetary + heater_power(i) - Q_radiation - Q_convection;
    
    % Temperature Change
    delta_T = Q_net * time_step / (payload_mass * payload_specific_heat);
    temperature_payload(i) = temperature_payload(i-1) + delta_T;
    
    % Clamp temperature
    temperature_payload(i) = max(min(temperature_payload(i), target_temp_max), target_temp_min);
end

%% Plotting
figure('Position', [100, 100, 1200, 800]);

% Altitude Profile
subplot(3,1,1);
plot(time/60, altitude, 'b', 'LineWidth', 2);
title('Altitude Profile (Ascent and Descent)');
xlabel('Time (minutes)');
ylabel('Altitude (m)');
grid on;

% Temperature Evolution
subplot(3,1,2);
plot(time/60, temperature_payload, 'r', 'LineWidth', 2);
hold on;
plot(time/60, temperature_external, 'b--', 'LineWidth', 2);
yline(target_temp_min, 'g--', 'Min Temp');
yline(target_temp_max, 'm--', 'Max Temp');
title('Payload Temperature');
xlabel('Time (minutes)');
ylabel('Temperature (°C)');
legend('Payload Temp', 'External Temp', 'Min Temp', 'Max Temp');
grid on;

% Heater Power
subplot(3,1,3);
plot(time/60, heater_power, 'g', 'LineWidth', 2);
title('Heater Power');
xlabel('Time (minutes)');
ylabel('Power (W)');
grid on;

%% Results Summary
total_energy_consumed = trapz(time, heater_power);
fprintf('Total Heater Energy Consumed: %.2f J\n', total_energy_consumed);
fprintf('Maximum Heater Power: %.2f W\n', max(heater_power));
fprintf('Final Payload Temperature: %.2f °C\n', temperature_payload(end));
