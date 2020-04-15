
function [theta_out] = Simulate_InvertedPend_Tiff(theta_init, Tapp, simtime, dt, howmuch)

t = 0; %s
sim_length = simtime/dt;
% dt = dt*10e-3; %Convert milisecond dt to seconds

%Constants for simulating
% m = 84.3; %kg 
% h = .92; %m - center of mass height
% I = 85.2; %kg m^2
% g = 9.81;% m/s

I_factor = 0.243734;    % assumes form of Inertia calculation is m*h^I_factor
m_adjust = 0;    % percent change to mass
h_adjust = 0;       % percent change to height

% Peterka's base "normal subject" constants
m = 83.3; %kg (roughly 184 lb)
h = .896; %m - center of mass height (roughly 3 feet)
% I = 81.1; %kg m^2

m = m + m_adjust*m; %kg 
h = h + h_adjust*h; %m - center of mass height
I = m*h^I_factor; %kg m^2

g = 9.81;% m/s

theta = zeros(sim_length,2);
theta(1,:) = theta_init;    % in radians

for t = [1:1:sim_length]
    %Contribution from gravity
    Tg = m*g*h*sin(theta(t,1)); % state space (updating theta and theta dot simultaneously);
%     Tg = m*g*h*(theta(t,1));

    theta(t+1,:) = [theta(t,1) + theta(t,2)*dt + 1/2*(Tapp+Tg)/I*dt^2, theta(t,2) + (Tapp+Tg)/I*dt]; 
    %[position(initial P + initial V*dt + acceleration*dt^2), velocity(initial V + acceleration*dt)]

%     theta(t+1,1) = mod(theta(t+1,1)+pi,2*pi) - pi;  % wrap to pi
%     theta(t+1,1) = mod(theta(t+1,1)+10,20) - 10;  % wrap to 10 nA, equivalent to pi in simulation
end

if howmuch == 1
    theta_out = theta(end,:);
elseif howmuch == 2
    theta_out = theta(:,:);
end