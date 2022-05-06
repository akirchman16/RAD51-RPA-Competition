clearvars
close all;
tic
% This is the new model for the competition between RPA and RAD51
% nucleoproteins on ssDNA. It utilizes the Direct Gillespie Algorithm and
% tracks various saturations of the DNA lattice over time. RPA is split
% into two segments, A (which has a very strong affinity for ssDNA) and D
% (which has a lower affinity for ssDNA). Both proteins can bind and unbind
% whenever as long as there are available locations for it. Individual
% microscopic binding/unbinding is possible for both parts. Equilibrium is
% tested for using a linear fit to the last 1/4 of events. A linear fit
% with slope less than 0.01 (1% per unit time interval) and a y-int with
% less than 5% error with the average FracCover is considered the
% equilibrium qualification. No cooperativity is included between RPA
% molecules.

N = 5000;   %length of ssDNA lattice
% RAD51 Parameters
RAD51 = 51; %value that will be stored on lattice to represent bound RAD51
n_RAD51 = 3;    %length of RAD51 protein

L_RAD51_Total = 2;  %total concentration of RAD51 in solution
Percent_M_RAD51 = 0.5;    %percentage of RAD51 solution which is monomers
w_RAD51 = 1;    %cooperativity parameter for RAD51
k_on_RAD51 = 1;     %kinetic rate constant for RAD51 binding to ssDNA
k_off_RAD51 = 1;    %kinetic rate constant for RAD51 dissociating from ssDNA

L_RAD51_M = L_RAD51_Total*Percent_M_RAD51;  %calculates concentration of RAD51 monomers
L_RAD51_D = L_RAD51_Total-L_RAD51_M;    %calculates concentration of RAD51 dimers

% RPA Parameters
RPA_A = 1;  %value to represent A piece of RPA on lattice
RPA_D = 3;  %value to represent D piece of RPA on lattice
n_A = 10;   %length of A component of RPA
n_D = 10;   %length of D component of RPA

L_RPA = 2;  %concentration of RPA in solution
w_RPA = 1;  %cooperativity parameter of RPA (for macroscopic binding)
k_on_RPA_A = 50; %kinetic rate constant for RPA-A binding to ssDNA
k_on_RPA_D = 50;  %kinetic rate constant for RPA-D binding to ssDNA
k_off_RPA_A = 1; %kinetic rate constant for RPA-A dissociating from ssDNA
k_off_RPA_D = 1; %kinetic rate constant for RPA-D dissociating from ssDNA

n_RPA = sum([n_A,n_D]);   %calculates total length of RPA molecule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Memory Allocation
DNA = zeros(2,N);   %array to represent DNA
                    %top row is used to store locations of hinged open RPA-D
                    %bottom row actually represents the DNA itself
RAD51_Mon_BoundAtSpot = zeros(1,N); %array used to record where RAD51 Monomers are bound
RPA_A_BoundAtSpot = zeros(1,N); %array to record where RPA-A is actively bound
RPA_D_BoundAtSpot = zeros(1,N); %array to record where RPA-D is actively bound
RPA_D_HingedOpen = zeros(1,N);  %array to record where RPA-D is microscopically dissociated from lattice
RPA_A_HingedOpen = zeros(1,N);  %array to record where RPA-A is microscopically unbound
LocationHistory = zeros(13,1);  %Matrix used to store locations of all events. Same order as Full_Propensity

% Initial Values
t(1) = 0;   %initial time is zero
FracCover_RAD51(1) = 0; %initially no RAD51 is on the lattice
FracCover_RPA_A(1) = 0; %initially empty lattice
FracCover_RPA_D(1) = 0; %initially empty lattice
FracCover_RPA(1) = 0;   %initially no RPA on the lattice
FracCover_Total(1) = 0; %initially the lattice is empty
x_Bound_RAD51_M(1) = 0; %initially 0 RAD51 Monomers bound to lattice
x_Bound_RAD51_D(1) = 0; %initially 0 RAD51 Dimers bound to lattice
x_Bound_RPA_A(1) = 0;   %initially 0 RPA molecules bound to lattice
x_Bound_RPA_D(1) = 0;   %initially 0 RPA molecules bound to lattice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Event_Count = 0;    %counts how many events happen within the simulation
Equilibrium_RAD51 = 0;  %test of whether RAD51 saturation is at equilibrium (1 = at equilibrium)
Equilibrium_RPA = 0;    %test of whether RPA saturation is at equilibrium (1 = at equilibirium)
while any([Equilibrium_RAD51,Equilibrium_RPA] == 0) == 1 & t(end) <= 25  %runs the whole time that the system is not at equilibrium
    Event_Count = Event_Count+1;    %advances event counter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function: Search for Free Locations
    [FreeCounts,RPA_I,RPA_SC,RPA_DC,RAD51_Mon_I,RAD51_Mon_SC,RAD51_Mon_DC,RAD51_Dim_I,RAD51_Dim_SC,RAD51_Dim_DC] = LatticeSearch_Cluster(DNA,n_RAD51,n_A,n_D)
    %Bound Proteins/Hinged RPA-A and RPA-D
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Macro_Propensity = [k_on_RAD51*L_RAD51_M,k_on_RAD51*L_RAD51_M*w_RAD51,k_on_RAD51*L_RAD51_M*(w_RAD51^2),k_on_RAD51*L_RAD51_D,k_on_RAD51*L_RAD51_D*w_RAD51,k_on_RAD51*L_RAD51_D*(w_RAD51^2),k_on_RPA_A*L_RPA].*Populations;  %propensity functions for reactions
    Unbinding_Propensity = [k_off_RAD51*x_Bound_RAD51_M(Event_Count),k_off_RAD51*x_Bound_RAD51_D(Event_Count),k_off_RPA_A*x_Bound_RPA_A(Event_Count),k_off_RPA_D*x_Bound_RPA_D(Event_Count)];   %propensity functions of unbinding reactions
                % 1 - RAD51 Monomer Unbinding
                % 2 - RAD51 Dimer Unbinding
                % 3 - RPA-A Microscopically Dissociating
                % 4 - RPA-D Microscopically Dissociating
    Micro_Propensity = [k_on_RPA_A*numel(Free_for_RPA_A),k_on_RPA_D*numel(Free_for_RPA_D)];   %propensity function for RPA-D rebinding (hinging closed)
    Full_Propensity = [Macro_Propensity, Unbinding_Propensity, Micro_Propensity];  %full listing of propensity functions for the simulation
                % 1 - RAD51 Monomer I Binding
                % 2 - RAD51 Monomer SC Binding
                % 3 - RAD51 Monomer DC Binding
                % 4 - RAD51 Dimer I Binding
                % 5 - RAD51 Dimer SC Binding
                % 6 - RAD51 Dimer DC Binding
                % 7 - RPA Macro Binding
                % 8 - RAD51 Monomer Unbinding
                % 9 - RAD51 Dimer Unbinding
                % 10 - RPA-A Micro Dissociating
                % 11 - RPA-D Micro Dissociating
                % 12 - RPA-D Micro Binding
                % 13 - RPA-A Micro Binding
    
    a_0 = sum(Full_Propensity); %sum of all propensity functions
    Randoms = [rand,rand];  %random numbers for Monte Carlo step
    dt(Event_Count) = (1/a_0)*log(1/Randoms(1)); %time until next reaction occurs
    
    RPA = sort([RPA_I,RPA_SC,RPA_DC]);  %sorted list of all available RPA binding locations
    %Direct Method of Calculating which reaction is next to occur
    if Full_Propensity(1) > Randoms(2)*a_0     %RAD51 Monomer I Binding Occurs
        j(Event_Count) = 1;  %records which reaction occured
        RAD51_Mon_I_Bind_Spot = RAD51_Mon_I(randi(numel(RAD51_Mon_I)));    %selects random location to bind to that's available
        DNA(2,RAD51_Mon_I_Bind_Spot:RAD51_Mon_I_Bind_Spot+(n_RAD51-1)) = RAD51; %binds RAD51 monomer to DNA lattice
        RAD51_Mon_BoundAtSpot(RAD51_Mon_I_Bind_Spot) = 1;   %records that there is a protein bound at this location
        LocationHistory(1,Event_Count) = RAD51_Mon_I_Bind_Spot;   %stores location in LocationHistory
    elseif sum(Full_Propensity(1:2)) > Randoms(2)*a_0  %RAD51 Monomer SC Binding Occurs
        j(Event_Count) = 2; %records which reaction occured
        RAD51_Mon_SC_Bind_Spot = RAD51_Mon_SC(randi(numel(RAD51_Mon_SC)));  %selects random location to bind RAD51 to that's available for corresponding reaction
        DNA(2,RAD51_Mon_SC_Bind_Spot:RAD51_Mon_SC_Bind_Spot+(n_RAD51-1)) = RAD51;   %binds RAD51 monomer
        RAD51_Mon_BoundAtSpot(RAD51_Mon_SC_Bind_Spot) = 1;  %records that a RAD51 Mon. is bound here
        LocationHistory(2,Event_Count) = RAD51_Mon_SC_Bind_Spot;    %records location of the event
    elseif sum(Full_Propensity(1:3)) > Randoms(2)*a_0  %RAD51 Monomer DC Binding Occurs
        j(Event_Count) = 3; %records which reaction occured
        RAD51_Mon_DC_Bind_Spot = RAD51_Mon_DC(randi(numel(RAD51_Mon_DC)));  %selects random location for binding
        DNA(2,RAD51_Mon_DC_Bind_Spot:RAD51_Mon_DC_Bind_Spot+(n_RAD51-1)) = RAD51;   %binds RAD51 dimer to location
        RAD51_Mon_BoundAtSpot(RAD51_Mon_DC_Bind_Spot) = 1; %records where monomer of RAD51 is bound
        LocationHistory(3,Event_Count) = RAD51_Mon_DC_Bind_Spot;    %records where binding happened and when
    elseif sum(Full_Propensity(1:4)) > Randoms(2)*a_0  %RAD51 Dimer I Binding Occurs
        j(Event_Count) = 4; %records which reaction occured
        RAD51_Dim_I_Bind_Spot = RAD51_Dim_I(randi(numel(RAD51_Dim_I))); %selects random location for dimer binding
        DNA(2,RAD51_Dim_I_Bind_Spot:RAD51_Dim_I_Bind_Spot+(2*n_RAD51-1)) = RAD51;   %binds RAD51 dimer
        RAD51_Mon_BoundAtSpot([RAD51_Dim_I_Bind_Spot,RAD51_Dim_I_Bind_Spot+n_RAD51]) = 1;   %records location of each monomer in the dimer
        LocationHistory(4,Event_Count) = RAD51_Dim_I_Bind_Spot; %records where location happened with this event
    elseif sum(Full_Propensity(1:5)) > Randoms(2)*a_0  %RAD51 Dimer SC Binding Occurs
        j(Event_Count) = 5; %records which reaction occured
        RAD51_Dim_SC_Bind_Spot = RAD51_Dim_SC(randi(numel(RAD51_Dim_SC)));  %selects location for dimer binding
        DNA(2,RAD51_Dim_SC_Bind_Spot:RAD51_Dim_SC_Bind_Spot+(2*n_RAD51-1)) = RAD51; %binds dimer to location
        RAD51_Mon_BoundAtSpot([RAD51_Dim_SC_Bind_Spot,RAD51_Dim_SC_Bind_Spot+n_RAD51]) = 1; %reocrds location of each monomer part
        LocationHistory(5,Event_Count) = RAD51_Dim_SC_Bind_Spot;    %records where binding happened and which type of reaction
    elseif sum(Full_Propensity(1:6)) > Randoms(2)*a_0  %RAD51 Dimer DC Binding Occurs
        j(Event_Count) = 6; %records which reaction occured
        RAD51_Dim_DC_Bind_Spot = RAD51_Dim_DC(randi(numel(RAD51_Dim_DC)));  %selects location for dimer binding
        DNA(2,RAD51_Dim_DC_Bind_Spot:RAD51_Dim_DC_Bind_Spot+(2*n_RAD51-1)) = RAD51; %binds dimer to location
        RAD51_Mon_BoundAtSpot([RAD51_Dim_DC_Bind_Spot,RAD51_Dim_DC_Bind_Spot+n_RAD51]) = 1; %reocrds location of each monomer part
        LocationHistory(6,Event_Count) = RAD51_Dim_DC_Bind_Spot;    %records where binding happened and which type of reaction
    elseif sum(Full_Propensity(1:7)) > Randoms(2)*a_0  %RPA Macro I Binding Occurs
        j(Event_Count) = 7; %records which reaction occured
        RPA_Macro_Bind_Spot = RPA(randi(numel(RPA))); %chooses random location for RPA macro binding
        DNA(2,RPA_Macro_Bind_Spot:RPA_Macro_Bind_Spot+(n_A-1)) = RPA_A; %binds A part of RPA
        DNA(2,RPA_Macro_Bind_Spot+n_A:RPA_Macro_Bind_Spot+n_A+(n_D-1)) = RPA_D; %binds D component of RPA to lattice
        RPA_A_BoundAtSpot(RPA_Macro_Bind_Spot) = 1;   %records where A component is bound
        RPA_D_BoundAtSpot(RPA_Macro_Bind_Spot+n_A) = 1;  %records wehre D compoenet of RPA is bound
        LocationHistory(7,Event_Count) = RPA_Macro_Bind_Spot; %records where this reaction occured in this event
    elseif sum(Full_Propensity(1:8)) > Randoms(2)*a_0  %RAD51 Monomer Unbinding
        j(Event_Count) = 8; %records which reaction occured
        RAD51_Bound_Monomers = find(RAD51_Mon_BoundAtSpot == 1);    %list of all locations a RAD51 monomer is bound
        RAD51_M_Unbind_Spot = RAD51_Bound_Monomers(randi(numel(RAD51_Bound_Monomers))); %selects random monomer that will unbind
        DNA(2,RAD51_M_Unbind_Spot:RAD51_M_Unbind_Spot+(n_RAD51-1)) = 0; %unbinds the monomer
        RAD51_Mon_BoundAtSpot(RAD51_M_Unbind_Spot) = 0; %removes location from BountAtSpot
        LocationHistory(8,Event_Count) = RAD51_M_Unbind_Spot;  %records where this event occured
    elseif sum(Full_Propensity(1:9)) > Randoms(2)*a_0  %RAD51 Dimer Unbinding Occurs
        j(Event_Count) = 9; %RAD51 Dimer Unbinding
        RAD51_Bound_Dimers = Left_RAD51_Dimer_Filament; %list of all locations where a dimer is bound
        RAD51_D_Unbind_Spot = RAD51_Bound_Dimers(randi(numel(RAD51_Bound_Dimers))); %randomly selects dimer to unbind from lattice
        DNA(2,RAD51_D_Unbind_Spot:RAD51_D_Unbind_Spot+(2*n_RAD51-1)) = 0;   %unbinds the dimer
        RAD51_Dim_BoundAtSpot(RAD51_D_Unbind_Spot) = 0; %removes location from BoundAtSpot (not necessary)
        RAD51_Mon_BoundAtSpot([RAD51_D_Unbind_Spot,RAD51_D_Unbind_Spot+n_RAD51]) = 0; %removes location from RAD51 Mon. BoundAtSpot record array
        LocationHistory(9,Event_Count) = RAD51_D_Unbind_Spot;  %records where this event occured
    elseif sum(Full_Propensity(1:10)) > Randoms(2)*a_0 %RPA-A Micro Dissociation
        j(Event_Count) = 10; %records which reaction occured
%         RPA_A_CurrentBound = find(RPA_A_BoundAtSpot == 1);  %all locations where RPA-A is currently bound
        RPA_A_Unbind_Spot = A_Hingeable(randi(numel(A_Hingeable)));    %location of RPA_A unbinding
        DNA(2,RPA_A_Unbind_Spot:RPA_A_Unbind_Spot+(n_A-1)) = 0;  %hinges RPA-A open
        DNA(1,RPA_A_Unbind_Spot:RPA_A_Unbind_Spot+(n_A-1)) = RPA_A;
        RPA_A_BoundAtSpot(RPA_A_Unbind_Spot) = 0;   %removes location from BoundAtSpot
        LocationHistory(10,Event_Count) = RPA_A_Unbind_Spot;    %records where event occured
        if DNA(1,RPA_A_Unbind_Spot+n_A) == RPA_D   %if both RPA-A and RPA-D are hinged open...
            DNA(1,RPA_A_Unbind_Spot:RPA_A_Unbind_Spot+n_RPA-1) = 0; %...remove it from the DNA lattice
        end
    elseif sum(Full_Propensity(1:11)) > Randoms(2)*a_0 %RPA-D Hinging Open Occurs
        j(Event_Count) = 11; %records which reaction occured
%         RPA_D_CurrentBound = find(RPA_D_BoundAtSpot == 1);  %all locations where RPA-A is currently bound
        RPA_D_Unbind_Spot = D_Hingeable(randi(numel(D_Hingeable)));    %location of RPA_A unbinding
        DNA(2,RPA_D_Unbind_Spot:RPA_D_Unbind_Spot+(n_D-1)) = 0;  %hinges RPA-A open
        DNA(1,RPA_D_Unbind_Spot:RPA_D_Unbind_Spot+(n_D-1)) = RPA_D;
        RPA_D_BoundAtSpot(RPA_D_Unbind_Spot) = 0;   %removes location from BoundAtSpot
        LocationHistory(11,Event_Count) = RPA_D_Unbind_Spot;    %records where event occured
        if DNA(1,RPA_D_Unbind_Spot-1) == RPA_A   %if both RPA-A and RPA-D are hinged open...
            DNA(1,RPA_D_Unbind_Spot:RPA_D_Unbind_Spot+n_D-1) = 0; %...remove it from the DNA lattice
            DNA(1,RPA_D_Unbind_Spot-n_A:RPA_D_Unbind_Spot-1) = 0;
        end
    elseif sum(Full_Propensity(1:12)) > Randoms(2)*a_0 %RPA-D Micro Binding
        j(Event_Count) = 12; %records which reaction occured
        RPA_D_Micro_Bind_Spot = Free_for_RPA_D(randi(numel(Free_for_RPA_D)));  %random location where RPA-D can rebind (and will in this event)
        DNA(2,RPA_D_Micro_Bind_Spot:RPA_D_Micro_Bind_Spot+(n_D-1)) = RPA_D; %binds RPA-D to the DNA lattice
        DNA(1,RPA_D_Micro_Bind_Spot:RPA_D_Micro_Bind_Spot+(n_D-1)) = 0; %clears RPA-D from being hinged open
        RPA_D_BoundAtSpot(RPA_D_Micro_Bind_Spot) = 1;   %records that RPA-D is now bound to the DNA lattice
        RPA_D_HingedOpen(RPA_D_Micro_Bind_Spot) = 0;    %records that RPA-D is no longer hinged open
        LocationHistory(12,Event_Count) = RPA_D_Micro_Bind_Spot;    %records where this event occured
    elseif sum(Full_Propensity(1:13)) > Randoms(2)*a_0 %RPA-A Micro Binding
        j(Event_Count) = 13; %records which reaction occured
        RPA_A_Micro_Bind_Spot = Free_for_RPA_A(randi(numel(Free_for_RPA_D)));   %random location where RPA-D can rebind to the lattice
        DNA(2,RPA_A_Micro_Bind_Spot:RPA_A_Micro_Bind_Spot+(n_A-1)) = RPA_A; %binds RPA-A to the DNA lattice
        DNA(1,RPA_A_Micro_Bind_Spot:RPA_A_Micro_Bind_Spot+(n_A-1)) = 0; %removes RPA-A from the hinged open position
        RPA_A_BoundAtSpot(RPA_A_Micro_Bind_Spot) = 1;   %records that RPA-A is now bound to the lattice
        RPA_A_HingedOpen(RPA_A_Micro_Bind_Spot) = 0;    %records that RPA-A is no longer hinged open
        LocationHistory(13,Event_Count) = RPA_A_Micro_Bind_Spot;    %records where the RPA-A micro binding occured
    end
    
    %Bug Checking - Finding Whole Proteins
    if numel(find(DNA(2,:) == RAD51))/n_RAD51 ~= round(numel(find(DNA(2,:) == RAD51))/n_RAD51)  %checks if there are an integer number of RAD51 proteins
        disp('ERROR: BROKEN RAD51');
        disp(['Last Reaction: ', num2str(j(end))]);
        break;
    elseif numel(find(DNA(2,:) == RPA_A))/n_A ~= round(numel(find(DNA(2,:) == RPA_A))/n_A)     %checks if there is an integer number of RPA-A bound
        disp('ERROR: BROKEN RPA-A');
        disp(['Last Reaction: ', num2str(j(end))]);
        break;
    elseif numel(find(DNA(2,:) == RPA_D))/n_D ~= round(numel(find(DNA(2,:) == RPA_D))/n_D)      %checks for an integer number of bound RPA_D
        disp('ERROR: BROKEN RPA-D (Bound)');
        disp(['Last Reaction: ', num2str(j(end))]);
        break;
    elseif numel(find(DNA(1,:) == RPA_D))/n_D ~= round(numel(find(DNA(1,:) == RPA_D))/n_D)      %checks for an integer number of hinged open RPA-D
        disp('ERROR: BROKEN RPA-D (Open)');
        disp(['Last Reaction: ', num2str(j(end))]);
        break;
    end
        
    t(Event_Count+1) = t(Event_Count)+dt(Event_Count);  %advance time according to time interval selected by Gillespie Algorithm
    
    FracCover_RAD51(Event_Count+1) = numel(find(DNA(2,:) == RAD51))/N; %RAD51 saturation of the DNA lattice after each event
    FracCover_RPA_A(Event_Count+1) = numel(find(DNA(2,:) == RPA_A))/N;  %saturation of RPA-A
    FracCover_RPA_D(Event_Count+1) = numel(find(DNA(2,:) == RPA_D))/N;  %saturation of RPA-D on the lattice
    FracCover_RPA(Event_Count+1) = sum([FracCover_RPA_A(Event_Count+1),FracCover_RPA_D(Event_Count+1)]);  %saturation of all parts of RPA
    FracCover_Total(Event_Count+1) = sum([FracCover_RAD51(Event_Count+1),FracCover_RPA(Event_Count+1)]);  %total saturation of the lattice (both RPA and RAD51)
    
% Equilibrium Testing - Linear Slope & Intercept Method %%%%%%%%%%%%%%%%%%%
    if Event_Count >= 1000 %only tests for equilibrium after 1000 events have occured
        t_Equilibrium_Test = t(Event_Count+1-round(0.25*(Event_Count+1)):end);  %time values that we're testing for equilibrium
        RPA_Equilibrium_Test = FracCover_RPA(((Event_Count+1)-round(0.25*(Event_Count+1))):end);   %last 1/4 of Events saturation data for RPA
        RAD51_Equilibrium_Test = FracCover_RAD51(((Event_Count+1)-round(0.25*(Event_Count+1))):end);   %last 1/4 of Events saturation data for RAD51
        
        RPA_Avg_Saturation = sum(RPA_Equilibrium_Test)/numel(RPA_Equilibrium_Test); %average saturation in last 1/4 of Events (RPA)
        RAD51_Avg_Saturation = sum(RAD51_Equilibrium_Test)/numel(RAD51_Equilibrium_Test);   %average saturation in last 1/4 of Events (RAD51)
        
        RAD51_Fit = polyfit(t_Equilibrium_Test,RAD51_Equilibrium_Test,1);   %linear fit for RAD51 data (slope, y-int)
        RAD51_Yint_Error = abs(RAD51_Avg_Saturation-RAD51_Fit(2))/RAD51_Avg_Saturation; %y-intercept error of linear fit for RAD51 data
        RPA_Fit = polyfit(t_Equilibrium_Test,RPA_Equilibrium_Test,1);   %linear fit to RPA data (slope, y-int)
        RPA_Yint_Error = abs(RPA_Avg_Saturation-RPA_Fit(2))/RPA_Avg_Saturation;    %y-intercept error compared to average RPA saturation
        
        if abs(RPA_Fit(1)) < 0.01 & (RPA_Yint_Error < 0.05 | isnan(RPA_Yint_Error))   %if slope of RPA data is essentially zero and y-intercept is very close to avg. saturation value... (slope limit is change in saturation of 1% (~17 proteins) per 1 time interval)
            Equilibrium_RPA = 1;    %...then at equilibrium
        else
            Equilibrium_RPA = 0;    %...otherwise reset to not at equilibrium
        end
        if abs(RAD51_Fit(1)) < 0.01 & (RAD51_Yint_Error < 0.05 | isnan(RAD51_Yint_Error)) %if the slope of RAD51 data is essentially zero and y-intercept is very close to avg. saturation value... (slope limit is change in saturation of 1% (~3 proteins) per 1 time interval)
            Equilibrium_RAD51 = 1;  %...then we're at equilibrium
        else
            Equilibrium_RAD51 = 0;    %...otherwise reset to not at equilibrium
        end
    end
end

t_Equilibrium = t(Event_Count+1-round(0.25*(Event_Count+1)));   %time where equilibrium occured
toc
figure(1);  %plots of saturation over time
scatter(t,FracCover_RAD51,1,'red','filled');    %RAD51 Saturation
hold on;
scatter(t,FracCover_RPA_A,1,'cyan','filled');   %RPA-A Saturation
scatter(t,FracCover_RPA_D,1,'blue','filled');   %RPA-D Saturation
scatter(t,FracCover_RPA,1,'magenta','filled');  %RPA Saturation
scatter(t,FracCover_Total,1,'black','filled');  %total protein saturation
xline(t_Equilibrium, '--k',['t: ', num2str(round(t_Equilibrium,2))],'LabelHorizontalAlignment','left');
yline(RPA_Avg_Saturation,'--k',['RPA: ', num2str(round(RPA_Avg_Saturation,2))],'LabelHorizontalAlignment','left');
yline(RAD51_Avg_Saturation,'--k',['RAD51: ', num2str(round(RAD51_Avg_Saturation,2))],'LabelHorizontalAlignment','left');
xlabel('Time, t');
xlim([0 max(t)]);
ylabel('Saturation');
ylim([0 1]);
title('RAD51/RPA Competition Saturation');
legend('RAD51','RPA-A','RPA-D','All RPA','Total','location','southoutside','orientation','horizontal');
box on;

disp(['RAD51 Saturation: ', num2str(round(RAD51_Avg_Saturation,2))]);
disp(['RPA Saturation: ', num2str(round(RPA_Avg_Saturation,2))]);
disp(['Equilibrium Time: ', num2str(round(t_Equilibrium,2))]);