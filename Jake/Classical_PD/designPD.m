%this program design the PD controller for balance control using the
%inverted pendulum model
%Author: Jake Chung
%Spring 2020 (COVID19)

clear
close all
clc

%% Define physical parameters.
m = 82; %mass in kg
R = 0.05; %radius of the cylinder in meters
L = 1.83; %length of the cylinder in meters
phi = 80; %angle of the muscle attachment
a = 0.1; %length of the base of the cylinder to the attachment of the muscle m
b = 0.3; %damping of the inverted pendulum
I = 1/4*m*R^2+1/3*m*L^2;
d = a*sind(phi);
g = 9.81;

%% define the transfer function
s = tf('s');
Gp = 1/(I/d*s^2+b/d*s-m/d*g*L/2);

h = figure;
subplot(3,1,1)
margin(Gp)
grid on
%% Design the PD controller by finding the P and D constants
Kd = 17720;
Kp = 4.753;
KpKd_ratio = 2.061;
Gc = Kd*(KpKd_ratio + s);
G = feedback(Gc*Gp,1);

subplot(3,1,2)
margin(Gc)
grid on

subplot(3,1,3)
margin(Gc*Gp)
grid on

set(h,'Position',[1,41,700,1200]);

figure 
pzmap(Gp)

%% simulate the response
% endSim = 2;
% t = 0:0.001:endSim;
% u = zeros(1,length(t));
% u(1:round(1/5*length(t))) = 100;
% figure
% lsim(G,u,t)