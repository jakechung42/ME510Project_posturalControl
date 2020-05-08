%this program design the PD controller for balance control using the
%inverted pendulum model
%Author: Jake Chung
%Spring 2020 (COVID19)

clear
close all
clc

%% Define physical parameters.
m = 82; %mass in kg
R = 0.1; %radius of the cylinder in meters
L = 1.83; %length of the cylinder in meters
phi = 70; %angle of the muscle attachment
a = 0.1; %length of the base of the cylinder to the attachment of the muscle m
b = 0.3; %damping of the inverted pendulum
I = 1/4*m*R^2+1/3*m*L^2;
d = a*sind(phi);

%% define the transfer function
s = tf('s');
sys_tf = 1/(-I/d*s^2+b/d*s);

step(sys_tf)
