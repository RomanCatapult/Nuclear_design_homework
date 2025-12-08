%% Data definitions
clear;
close all;
clc;


% Geometrical data
barr_ext_diam                            = 2.5;                    % Barrel external diameter (m)
vess_int_diam                            = 3;                      % Vessel internal diameter (m)
ins_thick                                = 0.05;                   % Insulator thickness (m)
ins_thermal_cond                         = 1.4;                    % Thermal insulator conductivity (W/m/K)

% Primary fluid data
core_in_temp                             = 214 + 273.15;           % Core inlet temperature (K)
core_avg_out_temp                        = 254 + 273.15;           % Average core outlet temperature (K)
core_max_out_temp                        = 270 + 273.15;           % Maximum core outlet temperature (K)
p_in                                     = 7.5e6;                  % Internal pressure (Pa)
p_out                                    = 7.5e6;                  % External pressure (Pa)
mass_flow_rate                           = 3227;                   % Coolant mass flow rate (kg/s)
avg_specific_heat_prim                   = 4534;                   % Average specific heat coefficient for the primary fluid (J/kg/K)
avg_dens_prim                            = 852.5;                  % Average density for the primary fluid (kg/m^3)
avg_din_visc_prim                        = 1.259e-4;               % Average dinamic viscosity for the primary fluid (Pa*s)
avg_therm_conduc_prim                    = 0.658;                  % Average thermal conductivity for the primary fluid (W/m/K)

% Containment
cont_temp                                = 70 + 273.15;            % Containment water temperature
cont_press                               = 7.5e6;                  % Containment water pressure
avg_specific_heat_cont                   = 4172.5;                 % Average specific heat coefficient for the containment water (J/kg/K)
avg_dens_cont                            = 981.2;                  % Average density for the containment fluid (kg/m^3)
avg_din_visc_cont                        = 4.06e-4;                % Average dinamic viscosity for the containment fluid (Pa*s)
avg_therm_conduc_cont                    = 0.666;                  % Average thermal conductivity for the containment fluid (W/m/K)
thermal_exp_coeff_cont                   = 5.57e-4;                % Thermal expansion coefficient for the containment fluid (1/K)

% Steel
young                                    = 177e9;                  % Young's modulus (Pa)
poisson_coeff                            = 0.3;                    % poisson coefficient
linear_th_exp_coeff                      = 1.7e-5;                 % Linear thermal expansion coefficient (1/K)
th_cond_ves                              = 48.1;                   % Thermal conductivity (W/m/K)
mu_attenuation                           = 24;                     % Effective linear attenuation coefficient (1/m)

% Radiation source
flux_0                                   = 1.5e13;                 % Photon flux that impacts on the internal surface of the vessel (photons/cm^2/s)
E_ph                                     = 6 * 1e6 * 1.602e-19;    % Average photon energy (J/photon)
B                                        = 1.4;                    % Build-up factor

% Other data (not directly from table)
sigma_T_tab                              = 0.86;
p_design                                 = p_in * 1.1;             % Design pressure with 10 % safety margin (Pa)
T_design_vess                            = 350;                    % Design temperature, from which we get Sm (°C)
core_height                              = 4;                      % Core height (not vessel height!)
kinematic_viscosity                      = (4.06e-4) / 981.2;      % Kinematic viscosity
Sm                                       = 114e6;                  % Value of Sm extracted considering T<350°C
Sy                                       = Sm * 3 / 2;             % Value for yield stress
T_melt                                   = 1450 + 273.15;          % Melting temperature of the steel is taken very low to simulate the worst case scenario(K)


%% Function definitions

% Adimensional numbers
Reynolds = @(rho, v, D, mu) (rho * v * D) / mu;                                 % Reynolds function
Prandtl = @(mu, Cp, k) (mu * Cp) / k;                                           % Prandtl function
Grashof = @(beta, DeltaT, L, nu) (9.81 * beta * DeltaT * (L^3)) / (nu^2);       % Grashof function
Dittus_boelter = @(Re, Pr) 0.023 * (Re^0.8) * (Pr^0.4);                         % Dittus-Boelter function for nusselt's number
Mc_Adams = @(Gr, Pr) 0.13 * ((Gr * Pr)^(1/3));                                  % McAdams function for nusselt's number
h = @(Nu, k, diameter) (Nu * k) / diameter;                                     % Function for convective heat exchange coefficient

% Temperature functions
C1 = @(T1, T2, q0, k_v, mu, L, h1, u2) ((q0 / (k_v * mu^2)) * (exp(-mu * L)-1) - (q0 / mu) * (1 / h1 + exp(-mu * L) / u2) - (T1 - T2)) / (L + k_v/h1 + k_v/u2);
C2 = @(T1, q0, h1, mu, C1, k_v) T1 + q0 / (h1 * mu) + C1 * k_v / h1 + q0 / (k_v * mu^2);
T_fun = @(T1, T2, L, h1, u2, x, qv, mu, k_v) -qv / (k_v * mu^2) .* exp(-mu .* x) + C1(T1, T2, qv, k_v, mu, L, h1, u2) .* x + C2(T1, qv, h1, mu, C1(T1, T2, qv, k_v, mu, L, h1, u2), k_v);
max_thermal_stress_fun = @(sigma_T, alpha, E, q_V, th_cond, nu, mu_attenuation) (sigma_T * alpha * E * q_V) / (th_cond * (1 - nu) * mu_attenuation^2);

% Stress functions
sigma_r_mar_fun = @(p) - p / 2;
sigma_theta_mar_fun = @(p, d, s) (p * (d / 2)) / s;
sigma_z_mar_fun = @(p, d, s) (p * (d / 2)) / (2 * s);
sigma_r_lame_fun = @(a, b, p_in, p_out, r) -((a^2 * b^2) / (b^2 - a^2)) * ((p_in - p_out) / (r^2)) + ((a^2 * p_in - b^2 * p_out) / (b^2 - a^2));
sigma_theta_lame_fun = @(a, b, p_in, p_out, r) ((a^2 * p_in - b^2 * p_out) / (b^2 - a^2)) + ((a^2 * b^2 * (p_in - p_out)) / (r^2 * (b^2 - a^2)));
sigma_z_lame_fun = @(sigma_r, sigma_theta, poisson_coeff) poisson_coeff * (sigma_r + sigma_theta);

% Corradi's pressures for buckling and plastic collapse verification (real tubes)
qE_corr_fun = @(young, poisson_coeff, D, t) ((2 * young) / (1 - poisson_coeff^2)) * 1 / ((D / t) * (D / t - 1)^2);
q0_corr_fun = @(Sy, D, t) 2 * Sy * (t / D) * (1 + t / (2 * D));


%% Vessel's thickness under design conditions
thickness_des_ves = (((p_design) * (vess_int_diam / 2)) / (Sm - (0.5 * p_design))) * 1.1;   % minimum thickness with 10 % safety margin. By multiplying by 1.1 i'm being double conservative, since i am already using p_design = p_in * 1.1 (m)


%% Buckling / plastic collapse verification

W = min(((vess_int_diam * 1000 + 50) / 200) * ((vess_int_diam * 1000 + 1250) / 200), vess_int_diam * 1000 / 100) / (vess_int_diam * 1000);                      % Maximum allowed ovality
t_corr = thickness_des_ves;    % Iterative thickness used for the corradi method, starting from the design thickness
s_corr = 2;
iter_corr = 0;

while true
    qE_corr = qE_corr_fun(young, poisson_coeff, vess_int_diam + t_corr, t_corr);
    q0_corr = q0_corr_fun(Sy, vess_int_diam + t_corr, t_corr);
    Z_corr = (sqrt(3) / 4) * (((2 * (vess_int_diam + t_corr)) / t_corr) + 1) * W;
    qU_corr = q0_corr / sqrt(1 + Z_corr^2);
    qL_corr = 0.5 * (q0_corr + qE_corr * (1 + Z_corr) - sqrt(((q0_corr + qE_corr * (1 + Z_corr))^2) - 4 * q0_corr * qE_corr));

    if (q0_corr / qE_corr) < 0.04
        mu_corr = 1;

    elseif (q0_corr / qE_corr) > 0.7
        mu_corr = 0;

    else
        mu_corr = 0.35 * log(qE_corr / q0_corr) - 0.125;

    end

    qC_corr = mu_corr * qU_corr + (1 - mu_corr) * qL_corr;
    p_all_corr = (1 / s_corr) * qC_corr;

    t_corr = t_corr + 0.001;                % Adding 1mm at each iteration
    iter_corr = iter_corr + 1;

    if p_all_corr > p_out
        break;
    end
end


%% Forced convection inside the vessel

A_ann = pi/4 * (vess_int_diam^2 - barr_ext_diam^2);
v_1  = mass_flow_rate / (avg_dens_prim * A_ann);

Re_1 = Reynolds(avg_dens_prim, v_1, 0.5, avg_din_visc_prim);
Pr_1 = Prandtl(avg_din_visc_prim, avg_specific_heat_prim, avg_therm_conduc_prim);


%% Convection inside and outside the vessel

% Case 1: between the inner surface of the vessel and the primary fluid
Nu_1 = Dittus_boelter(Re_1, Pr_1);                                                                  % Nusselt number
d_h1 = 0.5;                                                                                         % Hydraulic diameter
h1 = h(Nu_1, avg_therm_conduc_prim, d_h1);                                                          % Convective heat exchange coefficient

% Case 2: between the outer thermal insulation surface and the still water inside the containment (natural convection)
DeltaT_Gr_2 = 30;                                                                                   % Text says to assume this difference for Grashof number
Gr_2 = Grashof(thermal_exp_coeff_cont, DeltaT_Gr_2, core_height, kinematic_viscosity);              % Grashof number
Pr_2 = (avg_din_visc_cont * avg_specific_heat_cont) / avg_therm_conduc_cont;                        % Prandtl number
Nu_2 = Mc_Adams(Gr_2, Pr_2);                                                                        % Nusselt number
h2 = h(Nu_2, avg_therm_conduc_cont, core_height);                                                   % Convective heat exchange coefficient

% On the vessel outer surface
R_v_i = vess_int_diam / 2;                                                                          % Inner radius of the vessel
R_v_o = vess_int_diam / 2 + t_corr;                                                                 % Outer radius of the vessel
R_t_o = R_v_o + ins_thick;                                                                          % Outer radius up to the thermal insulator
thermal_resistence = ((R_v_o * log(R_t_o / R_v_o)) / ins_thermal_cond) + (R_v_o / (h2 * R_t_o));    % Thermal resistance between vessel and water in containment
u2 = 1 / thermal_resistence;                                                                        % Global heat exchange coefficient on the vessel outer surface (between vessel and water in containment)


%% Volumetric heat generation

I0 = B * flux_0 * E_ph * 1e4;                           % Energy flux on the internal surface of the vessel caused by gamma radiation (W/m^2)                            
q_V = zeros(1, 200);
x = linspace(0, t_corr, 200);

for i = 1:200
    q_V(i) = I0 * mu_attenuation * exp(- mu_attenuation * x(i));       % Internal heat generation in every point through the thickness (W/m^3)
end

figure;
plot(x, q_V/1e6, 'r');
grid on;
title("Radial profile of q''' across the vessel without thermal shield");
xlabel("x [m]");
ylabel("q''' [MW/m^3]");


%% Analytical temperature solution & thermal power flux at inner and outer surface of vessel with and without gamma radiation

x = linspace(0, t_corr, 200);
T_ves = zeros(1, 200);

for i = 1:200
    T_ves(i) = T_fun(core_in_temp, cont_temp, t_corr, h1, u2, x(i), q_V(1), mu_attenuation, th_cond_ves); % Temperature across the thickness of the vessel considering also gamma radiation (K)
end

T_ves = T_ves - 273.15;                                         % Temperature across the thickness of the vessel with gamma radiation (°C)
T_max_ves = max(T_ves);                                         % Maximum temperature across the vessel (°C)
x_T_max_ves = x(T_ves == max(T_ves));                           % Position in the thickness in which T is the highest (m from internal surface)     %%%%%%%%%% a davide viene leggermente inferiore
T_avg_vess = (1 / t_corr) * trapz(x, T_ves);                    % Average temperature in the vessel (°C)
T_avg_vess_vett = zeros(1, 200) + T_avg_vess;                   % Used to graph the average temperature in the vessel (°C)
creep_coeff_vess_design = (max(T_ves) + 273.15) / T_melt;       % If greater than 0.3 we need a thermal shield (in this case we do need it!)

% The design temperature we got is roughly 20 degrees more than the average
% temperature in the vessel, which is good considering that we still have to take into account the thermal shield 

qV_no_gamma = 0;                                                      % New value of power produced by gamma radiation, supposed null
T_ves_no_gamma = zeros(1, 200);

for i = 1:200
    T_ves_no_gamma(i) = T_fun(core_in_temp, cont_temp, t_corr, h1, u2, x(i), qV_no_gamma, mu_attenuation, th_cond_ves); 
end

T_ves_no_gamma = T_ves_no_gamma - 273.15;                               % Temperature across the thickness of the vessel without gamma radiation (°C)
q0_no_gamma = th_cond_ves * (T_ves_no_gamma(1) - T_ves_no_gamma(end));  % Heat flux without gamma generation (W/m)

figure;
plot(x, T_ves);
hold on;
plot(x, T_ves_no_gamma);
plot(x, T_avg_vess_vett)
grid on;
title("Temperature profile across the vessel thickness without thermal shield");
legend("T_\gamma", "T_{no \gamma}", "T_{media}");
xlabel("x [m]");
ylabel("T [°C]");


%% Mechanical stresses (Mariotte)

% First, using Mariotte's approach, we find the mechanical stresses across the vessel thickness
sigma_r_mar = zeros(1, 200);
sigma_theta_mar = zeros(1, 200);
sigma_z_mar = zeros(1, 200);
r = linspace(vess_int_diam / 2, vess_int_diam / 2 + t_corr, 200);

for i = 1:200
    sigma_r_mar(i) = sigma_r_mar_fun (p_design) / 1e6;                                                % Mariotte's radial stress component (MPa)
    sigma_theta_mar(i) = sigma_theta_mar_fun (p_design, vess_int_diam, t_corr) / 1e6;      % Mariotte's hoop stress component (MPa)
    sigma_z_mar(i) = sigma_z_mar_fun (p_design, vess_int_diam, t_corr) / 1e6;              % Mariotte's axial stress component (MPa)
end

sigma_mar_vett = [sigma_r_mar(1), sigma_theta_mar(1), sigma_z_mar(1)];                                % Vector containing mariotte's stresses (MPa)
sigma_c_prim_tresca = max([abs(sigma_mar_vett(1) - sigma_mar_vett(2)), ...
    abs(sigma_mar_vett(2) - sigma_mar_vett(3)), abs(sigma_mar_vett(1) - sigma_mar_vett(3))]);         % Sigma_comparison found using Tresca's approach with Mariotte's stresses (MPa)

% sigma_c_prim_tresca should pass a preliminary test regarding primary stresses, which is that sigma_c_mar should be inferior to Sm, which it is in our case, so we are good


%% Mechanical stresses (Lame)

% Now, using Lamé's approach, we find the mechanical stresses across the vessel thickness to which we add the maximum mechanical stress
sigma_r_lame = zeros(1, 200);
sigma_theta_lame = zeros(1, 200);
sigma_z_lame = zeros(1, 200);
p_out_lame = 0;                         % If using the same value as p_in, the we would get sigma_r_M_lame = sigma_theta_M_lame (Pa)

for i = 1:200
    sigma_r_lame(i) = sigma_r_lame_fun (vess_int_diam / 2, (vess_int_diam / 2) + t_corr, p_design, p_out_lame, r(i)) / 1e6;              % Lamé's radial stress component (MPa)
    sigma_theta_lame(i) = sigma_theta_lame_fun (vess_int_diam / 2, (vess_int_diam / 2) + t_corr, p_design, p_out_lame, r(i)) / 1e6;      % Lamé's hoop stress component (MPa)
    sigma_z_lame(i) = sigma_z_lame_fun (sigma_r_lame(i), sigma_theta_lame(i), poisson_coeff);                                                       % Lamé's axial stress component (MPa)
end

figure;
plot(x, sigma_r_mar);
hold on;
plot(x, sigma_r_lame);
plot(x, sigma_theta_mar);
plot(x, sigma_theta_lame);
plot(x, sigma_z_mar);
plot(x, sigma_z_lame);
grid on;
legend("\sigma_r^{mar}", "\sigma_r^{lame}", "\sigma_\theta^{mar}", "\sigma_\theta^{lame}", "\sigma_z^{mar}", "\sigma_z^{lame}");
title("Mariotte's and Lamé's stresses across the vessel's thickness");
xlabel("x [m]");
ylabel("Stress [MPa]")

max_thermal_stress = max_thermal_stress_fun (sigma_T_tab, linear_th_exp_coeff, young, q_V(1), th_cond_ves, poisson_coeff, mu_attenuation) / 1e6;       % Maximum thermal stress (MPa)
sigma_lame_vett = [sigma_r_lame(1), sigma_theta_lame(1), sigma_z_lame(1)];                                                                             %% qui sto usando il vettore sotto anzichè quello a questa riga per trovare sigma_c degli stress secondary
sigma_lame_thermal_vett = [sigma_r_lame(1), sigma_theta_lame(1) + max_thermal_stress, sigma_z_lame(1) + max_thermal_stress];                           % Sigma_comparison found using Tresca's approach with Lamé's stresses (MPa) 


%% Maximum stress

sigma_c_second_tresca = max([abs(sigma_lame_thermal_vett(1) - sigma_lame_thermal_vett(2)), abs(sigma_lame_thermal_vett(1) - sigma_lame_thermal_vett(3)), abs(sigma_lame_thermal_vett(2) - sigma_lame_thermal_vett(3))]); % Sigma_comparison using Tresca's approach with Lamé's approximation (MPa)

% If sigma_c_lame_tresca comes out to be greater than 3 Sm (which it does) then it means that we need a thermal shield!


%% Actual thermal stress across the vessels's thickness: radial and hoop components

sigma_theta_thermal = zeros(1, 200);
sigma_r_thermal = zeros(1, 200);

for i = 1:200
    sigma_theta_thermal(i) = ((linear_th_exp_coeff * young) / (1 - poisson_coeff)) * ((T_avg_vess - T_ves(i))) / 1e6;
end

figure;
plot(x, sigma_theta_thermal);
grid on;
legend("\sigma_\theta^{thermal}");
title("\sigma_\theta^{thermal} profile across the vessel's thickness");
xlabel("x [m]");
ylabel("Thermal Stress [MPa]")

sigma_theta_thermal_max = max(sigma_theta_thermal);
sigma_z_thermal_max = sigma_theta_thermal_max;

% In order to be conservative, we considered the maximum thermal stress
% found using the formula given by the homework file


%% Thermal shield

% In order to calculate the thickness of the thermal shield, I must first
% evaluate the maximum allowable secondary stress in the vessel, and from
% that I must reverse engineer the thickness needed for the thermal shield.
% The maximum allowable secondary stress must be below Sm and in order to
% be conservative, I shall take 30 % less than that.

max_all_second_stress_vess = 0.7 * 3 * Sm / 1e6;

% The new value for sigma_c_second_tresca must be lower
% max_all_second_stress_vess, therefore i can start an iterative cycle in
% which I increase the thermal shield's thickness and stop when the new
% value for sigma_c_second_tresca is lower than max_all_second_stress_vess.

% Also, the thickness of the thermal shield shall lower the temperature of
% the vessel enough to exit the thermal creep's region, so we must chech
% the average temperature of the vessel ad each iteration and verify that
% it guarantees no thermal creep

iter = 0;
thickness_th_sh = 0;
T_ves_with_th_shield = zeros(1, 200);
T_th_shield = zeros(1, 200);

while true
    iter = iter + 1;
    thickness_th_sh = thickness_th_sh + 0.001; % Increasing the thermal shield's thickness by 1mm each iteration
    x_th_sh = linspace(0, thickness_th_sh, 200);
    q0_V_vess_with_th_sh = I0 * mu_attenuation * exp(- mu_attenuation * thickness_th_sh);
    max_thermal_stress_with_th_sh = max_thermal_stress_fun(sigma_T_tab, linear_th_exp_coeff, young, q0_V_vess_with_th_sh, th_cond_ves, poisson_coeff, mu_attenuation) / 1e6;
    sigma_lame_thermal_vett_with_th_shield = [sigma_r_lame(1), sigma_theta_lame(1) + max_thermal_stress_with_th_sh, sigma_z_lame(1) + max_thermal_stress_with_th_sh];
    sigma_c_second_tresca_with_th_shield = max([abs(sigma_lame_thermal_vett_with_th_shield(1) - sigma_lame_thermal_vett_with_th_shield(2)), abs(sigma_lame_thermal_vett_with_th_shield(1) - sigma_lame_thermal_vett_with_th_shield(3)), abs(sigma_lame_thermal_vett_with_th_shield(2) - sigma_lame_thermal_vett_with_th_shield(3))]);
    

    % Geometrical data
    sh_center = (vess_int_diam - (vess_int_diam - barr_ext_diam) / 2) / 2;
    R_s_i = sh_center - thickness_th_sh / 2;
    R_s_o = sh_center + thickness_th_sh / 2;


    % Annulus between barrel and internal surface of thermal shield
    A1_sh = pi * (R_s_i^2 - (barr_ext_diam / 2)^2);
    v1_sh = mass_flow_rate / (avg_dens_prim * A1_sh);

    Re_1_sh = Reynolds(avg_dens_prim, v1_sh, 2 * R_s_i - 2.5, avg_din_visc_prim);
    Pr_1_sh = Prandtl(avg_din_visc_prim, avg_specific_heat_prim, avg_therm_conduc_prim);
    Nu_1_sh = Dittus_boelter(Re_1_sh, Pr_1_sh);
    h1_sh = h(Nu_1_sh, avg_therm_conduc_prim, 2 * R_s_i - 2.5);


    % Annulus between external thermal shield and internal surface of vessel
    A2_sh = pi * ((vess_int_diam / 2)^2 - R_s_o^2);
    v2_sh = mass_flow_rate / (avg_dens_prim * A1_sh);

    Re_2_sh = Reynolds(avg_dens_prim, v2_sh, vess_int_diam / 2 - 2 * R_s_o, avg_din_visc_prim);
    Pr_2_sh = Prandtl(avg_din_visc_prim, avg_specific_heat_prim, avg_therm_conduc_prim);
    Nu_2_sh = Dittus_boelter(Re_2_sh, Pr_2_sh);
    h2_sh = h(Nu_2_sh, avg_therm_conduc_prim, 2 * vess_int_diam / 2 - 2 * R_s_o);

    % Temperature across the vessel's thickness with thermal shield
    for i = 1:200
        T_ves_with_th_shield(i) = T_fun(core_in_temp, cont_temp, t_corr, h1, u2, x(i), q0_V_vess_with_th_sh, mu_attenuation, th_cond_ves) - 273.15; 
    end
    T_avg_vess_with_th_shield = (1 / t_corr) * trapz(x, T_ves_with_th_shield);       % Average temperature in the vessel with thermal shield (°C)
    T_avg_vess_with_th_shield_vett = zeros(1, 200) + T_avg_vess_with_th_shield;

    % Temperature across the thermal shield's thickness
    for i = 1:200
        T_th_shield(i) = T_fun(core_in_temp, core_in_temp, thickness_th_sh, h1_sh, h1_sh, x_th_sh(i), (I0 * mu_attenuation), mu_attenuation, th_cond_ves) - 273.15; 
    end
    T_avg_th_shield = (1 / thickness_th_sh) * trapz(x_th_sh, T_th_shield);                      % Average temperature in the thermal shield (°C)
    T_avg_th_shield_vett = zeros(1, 200) + T_avg_th_shield;
    
    if sigma_c_second_tresca_with_th_shield < max_all_second_stress_vess && (max(T_ves_with_th_shield) + 273.15) < 0.29 * T_melt        % Conditions that arrests the cycle
        break;
    end
end

creep_coeff_with_th_shield = (max(T_ves_with_th_shield) + 273.15) / T_melt;
creep_coeff_th_shield = (max(T_th_shield) + 273.15) / T_melt;

% Plotting the temperature of the vessel with the thermal shield
figure;
plot(x, T_ves_with_th_shield);
hold on;
plot(x, T_avg_vess_with_th_shield_vett);
grid on;
legend("T_{ves}", "T_{avg, ves}");
title("Temperature inside the vessel with thermal shield");
xlabel("x [m]");
ylabel("T [°C]")

% Plotting the temperature of the thermal shield
figure;
plot(x_th_sh, T_th_shield);
hold on;
plot(x_th_sh, T_avg_th_shield_vett);
grid on;
legend("T_{th shield}", "T_{avg, th shield}");
title("Temperature inside the thermal shield");
xlabel("x [m]");
ylabel("T [°C]")


% Pensare al liner interno per la corrosione (forse)

% Pensare a cosa cambia se posizioniamo il thermal shield in diverse posizioni

% 