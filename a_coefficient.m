function [CA, CNa, Xcp] = a_coefficient(Mach)
    % Define Mach number data points
    M = [0 0.3 0.5 0.7 0.9 1.0 1.1 1.2 1.3 1.4 1.5 2.5 3.5 4.5 5.5];
    
    % Define corresponding aerodynamic coefficients and center of pressure
    xcp_row = [10.167 10.167 10.242 10.323 10.510 10.696 11.385 11.492 10.825 9.830 9.534 7.927 7.076 6.459 6.080];
    CA_row = [0.256 0.256 0.244 0.236 0.342 0.473 0.535 0.523 0.480 0.451 0.432 0.346 0.283 0.243 0.215];
    CNa_row = [0.0812 0.0812 0.0819 0.0826 0.0842 0.0872 0.0951 0.0962 0.1046 0.0968 0.0962 0.0831 0.0694 0.0606 0.0541];
    
    % Interpolate values based on the given Mach number
    Xcp = interp1(M, xcp_row, Mach, 'spline');
    CA = interp1(M, CA_row, Mach, 'spline');
    CNa = interp1(M, CNa_row, Mach, 'spline');
end
