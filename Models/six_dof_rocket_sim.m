function [dx, Jyy, Jzz, Fthrust, length_g] = six_dof_rocket_sim(t, x, u, P, a, rho, h, gt, gn, disturbance_p, disturbance_y)
    % Rocket parameters
    params = params_test();
    m_ignition = params.m_ignition; % Initial rocket mass
    m_cutoff = params.m_cutoff; % Mass at cutoff
    t_burn = params.t_burn; % Burn time
    Isp_sea_level = params.Isp_sea_level; % Specific impulse at sea level
    Isp_vacuum = params.Isp_vacuum; % Specific impulse in vacuum
    g0 = params.g0;
    Pressure_sea_level = params.Pressure_sea_level;
    Pressure_vaccume = params.Pressure_vaccume;
    Reference_length_m = params.Reference_length_m;
    Reference_area_m2 = params.Reference_area_m2;

    %% State variables
    Vx = x(1);
    Vy = x(2);
    Vz = x(3);
    u_b = x(4);
    v_b = x(5);
    w_b = x(6);
    phi = x(7);
    theta = 0.1;
    psi = x(9);
    omega = x(10:12);

    %% Control inputs
    deltaY = u(1);
    deltaP = u(2);
    %% Ensure values are finite and reasonable
    if any(isnan(x)) || any(isinf(x))
        error('State vector x contains NaN or Inf values');
    elseif any(isnan(u)) || any(isinf(u))
        error('Control input u contains NaN or Inf values');
    end

    %% DCM
    % Transformation to body frame
    R = euler_dcm(phi, theta, psi);

    %% Rocket Velocity
    Vwind = [0; 0; 0]; % Define wind velocity if it's supposed to be zero
    Vair = [Vx; Vy; Vz] - Vwind;
    Vair_norm = norm(Vair);
    Vbody = [u_b; v_b; w_b];
    Vbody_norm = norm(Vbody);

    %% Aerodynamics Parameters
    if Vair_norm <= 1e-10
        CA = 0;
        CN = 0;
        CS = 0;
        Xcp = 10.167;
        alpha = 0;
        beta = 0;
    else
        % Calculate Mach number
        Mach = Vbody_norm / a;

        % Mach vs. Aerodynamic Coefficient
        [CA, CNa, Xcp] = aero_coefficient(Mach);

        % Calculate AoA & flight-path angle
        alpha = atan(w_b / u_b);
        beta = asin(v_b / Vbody_norm);

        % Radians to Degrees
        alpha_d = rad2deg(alpha);
        beta_d = rad2deg(beta);

        % Calculate normal force coefficient CN and lateral slip force coefficient CS
        CN = CNa * alpha_d;
        CS = CNa * beta_d;
    end

    %% Aerodynamics Forces
    % Axial force (A)
    A = 0.5 * rho * Vair_norm^2 * Reference_area_m2 * CA;
    % Normal force (Na) dependent on AoA
    Na = 0.5 * rho * Vair_norm^2 * Reference_area_m2 * CN;
    % Lateral slip force (Sb) dependent on angle of slip
    Sb = 0.5 * rho * Vair_norm^2 * Reference_area_m2 * CS;
    % Consider Normal force & Lateral slip force
    N = Na * sin(alpha);
    S = Sb * sin(beta);
    % Resultant aerodynamic forces in the body frame
    Faero_body = [-A; S; -N];
    Faero = R * Faero_body;

    %% Thrust
    m_dot_avg = (m_ignition - m_cutoff) / t_burn; % Constant mass flow rate
    Isp = Isp_sea_level + ((Isp_vacuum - Isp_sea_level) / (Pressure_vaccume - Pressure_sea_level)) * (P - Pressure_sea_level);
    if t <= t_burn
        Fthrust = Isp * g0 * m_dot_avg;
    else
        Fthrust = 0;
    end
    Ft_body = Fthrust * [cos(deltaY) * cos(deltaP); -sin(deltaY); -cos(deltaY) * sin(deltaP)];
    Ft = R * Ft_body;

    %% Mass & Moment of Inertia
    if t == 0
        mass_kg = m_ignition;
        XCG = 8.37e-12 * m_ignition^3 - 2.80e-07 * m_ignition^2 + 3.05e-03 * m_ignition + 2.49;
    elseif t > 0 && t <= t_burn
        mass_kg = m_ignition - m_dot_avg * t;
        XCG = 8.37e-12 * mass_kg^3 - 2.80e-07 * mass_kg^2 + 3.05e-03 * mass_kg + 2.49;
    else
        mass_kg = m_cutoff;
        XCG = 8.37e-12 * m_cutoff^3 - 2.80e-07 * m_cutoff^2 + 3.05e-03 * m_cutoff + 2.49;
    end
    Jxx = -1.07e-06 * mass_kg^2 + 3.25e-02 * mass_kg + 4.72e+02;
    Jyy = 2.90e-07 * mass_kg^3 - 8.19e-03 * mass_kg^2 + 8.12e+01 * mass_kg - 2.62e+04;
    Jzz = 2.90e-07 * mass_kg^3 - 8.19e-03 * mass_kg^2 + 8.12e+01 * mass_kg - 2.62e+04;

    %% Calculate Lever Arm
    
    XCP = Xcp * Reference_length_m;
    length_g = -1 * abs(21 - XCG);
    length_a = -1 * abs(XCP - XCG);
    length_Aero = [length_a; 0; 0];
    length_GCM = [length_g; 0; 0];

    %% Gravity
    Fg = -1 * mass_kg * [gn; 0; gt];
    % Gravity (Body Frame)
    Fg_body = R' * Fg;

    %% Equation of Motion
    % Rotational motion
    Mt = -cross(Ft_body, length_GCM);
    Maero = -cross(Faero_body, length_Aero);
    M = Mt + Maero + [0; disturbance_p; disturbance_y];
    % Velocity equations
    delta_V = (Ft + Faero + Fg) / mass_kg ;
    delta_VB = ((Ft_body + Faero_body + Fg_body) / mass_kg ) - cross(omega, [u_b; v_b; w_b]);
    % Attitude equations
    delta_euler = [1, sin(phi) * tan(theta), cos(phi) * tan(theta);
        0, cos(phi), -sin(phi);
        0, sin(phi) / cos(theta), cos(phi) / cos(theta)] * omega;
    % Angular velocity equations
    delta_omega = zeros(3, 1);
    delta_omega(1) = 1 / Jxx * (M(1) - (Jyy - Jzz) * omega(2) * omega(3));
    delta_omega(2) = 1 / Jyy * (M(2) - (Jzz - Jxx) * omega(1) * omega(3));
    delta_omega(3) = 1 / Jzz * (M(3) - (Jxx - Jyy) * omega(1) * omega(2));

    % Check for non-finite values
    if any(isnan(delta_V)) || any(isinf(delta_V))
        error('Non-finite delta_V');
    elseif any(isnan(delta_euler)) || any(isinf(delta_euler))
        error('Non-finite delta_euler');
    elseif any(isnan(delta_VB)) || any(isinf(delta_VB))
        error('Non-finite delta_VB');
    elseif any(isnan(delta_omega)) || any(isinf(delta_omega))
        error('Non-finite delta_omega');
    end

    if h <= -1e-10 % Crashing condition
        dx = zeros(12,1);
        warning('Rocket Crashed !');
    else
        dx = [delta_V(1); delta_V(2); delta_V(3); delta_VB(1); delta_VB(2); delta_VB(3); delta_euler(1); delta_euler(2); delta_euler(3); delta_omega];
    end
end
