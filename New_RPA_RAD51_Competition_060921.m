clearvars
close all;

% This is the new model for the competition between RPA and RAD51
% nucleoproteins on ssDNA. It utilizes the Direct Gillespie Algorithm and
% tracks various saturations of the DNA lattice over time. RPA is split
% into two segments, A(which has a very strong affinity for ssDNA) and D
% (which has a lower affinity for ssDNA). Both proteins can bind and unbind
% whenever as long as there are available locations for it.

N = 8660;   %length of ssDNA lattice
% RAD51 Parameters
n_RAD51 = 3;    %length of RAD51 protein
L_RAD51_Total = 2;  %total concentration of RAD51 in solution
Percent_M_RAD51 = 1;    %percentage of RAD51 solution which is monomers
w_RAD51 = 1;    %cooperativity parameter for RAD51
k_on_RAD51 = 1;     %kinetic rate constant for RAD51 binding to ssDNA
k_off_RAD51 = 1;    %kinetic rate constant for RAD51 dissociating from ssDNA

L_RAD51_M = L_RAD51_Total*Percent_M_RAD51;  %calculates concentration of RAD51 monomers
L_RAD51_D = L_RAD51_Total-L_RAD51_M;    %calculates concentration of RAD51 dimers

% RPA Parameters
n_A = 10;   %length of A component of RPA
n_D = 10;   %length of D component of RPA
L_RPA = 1;  %concentration of RPA in solution
w_RPA = 1;  %cooperativity parameter of RPA (for macroscopic binding)
k_on_RPAa = 10; %kinetic rate constant for RPA-A binding to ssDNA
k_on_RPAd = 8;  %kinetic rate constant for RPA-D binding to ssDNA
k_off_RPAa = 1; %kinetic rate constant for RPA-A dissociating from ssDNA
k_off_RPAd = 1; %kinetic rate constant for RPA-D dissociating from ssDNA

n_RPA = sum(n_A,n_D);   %calculates total length of RPA molecule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DNA = zeros(2,N);   %array to represent DNA
                    %top row is used to store locations of hinged open RPA-D
                    %bottom row actually represents the DNA itself
MaxTime = 1.5;  %maximum time the simulation runs to

t(1) = 0;   %initial time is zero

Event_Count = 0;    %counts how many events happen within the simulation
while max(t) <= MaxTime
    Event_Count = Event_Count+1;    %advances event counter
end
