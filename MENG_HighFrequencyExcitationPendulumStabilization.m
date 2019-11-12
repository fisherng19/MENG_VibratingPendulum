% MENG Research - High Frequency Excitation Pendulum Stabilization
% Project: Modeling the stabilization of an inverted pendulum with vibrations
% Date: 10/08/2019 - 11/12/2019

%------------- BEGIN CODE --------------

clear all
close all
clc

%%
% Define the state equations of the governing differential equation of the system
f = @( t,z, gamma, mu, a, omega) ...
    [ z(2);
    - mu*z(2) + omega^2*sin(z(1)) + a*sin( z(1) - gamma )*cos(t) ];

% Select parameter values
omegaRes = 7;                                     % Position resolution
aRes = 6;                                     % Velocity resolution
mu = 0.10;                                      % Viscuous damping coefficient
a = linspace(0,0.05,aRes);                    % Dimensionless forcing amplitude
omega = linspace(0, 0.15, omegaRes);              % Ratio of excitation frequency to the natural frequency
gamma = 0;                                      % Direction of excitation

z0 = [0;0];                                     % Initial position and velocity values
posBasin = 4;
velBasin = 4;
stabilityBasin = zeros(posBasin,velBasin);      % Matrix for basin of attraction
tf = 4000;                                      % Integrating time interval
tol = 0.1;

for posVary = 1:size(stabilityBasin,1);
    for velVary = 1:size(stabilityBasin,1);

        stabilityVal.y1 = zeros(omegaRes,aRes);
        stabilityVal.y2 = zeros(omegaRes,aRes);
        ttlOptions = zeros(length(mu),length(a));

        % Solves the differential equation and writes the data to an end-value
        % stability matrix
        z0 = 0.005 .* [posVary, velVary];
        for i = 1:length(omega)
            for j = 1:length(a)
            options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'Events', @(t,z) odeEvent(t,z,gamma, mu, a(j), omega(i)));
            sol = ode45( @(t,z) f(t,z,gamma, mu, a(j), omega(i)), [0, tf], z0, options);
            stabilityVal.y1(i,j) = sol.y(1,end);
            stabilityVal.y2(i,j) = sol.y(2,end);    
            end
        end

        % Plots stability output based on input parameters
        contourf(omega,a,stabilityVal.y1');
        xlabel('omega value');
        ylabel('a Value');


        stabilityBasin(posVary,velVary) = sum(abs( stabilityVal.y1(:)) <= tol );
    end
end


%-------------- END CODE ---------------