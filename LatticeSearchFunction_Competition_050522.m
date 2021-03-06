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
% This will also be used to test the new LatticeSearch functions.

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
%   Function Output Values:
%     Free Counts - number of positions on the lattice that are free for each type of binding of each protein type
%           1: RPA Isolated Sites
%           2: RPA Singly Contiguous Sites
%           3: RPA Doubly Contiguous Sites
%           4: RAD51 Monomer Isolated Sites
%           5: RAD51 Monomer Singly Congiguous Sites
%           6: RAD51 Monomer Doubly Contiguous Sites
%           7: RAD51 Dimer Isolated Sites
%           8: RAD51 Dimer Singly Contiguous Sites
%           9: RAD51 Dimer Doubly Contiguous Sites
%     Bound Counts - number of bound proteins of each type that are able to
%     unbind, respectively
%           1: RAD51 Monomer
%           2: RAD51 Dimer
%           3: RPA-A Hingeable
%           4: RPA-A Bindable
%           5: RPA-D Hingeable
%           6: RPA-D Bindable
   [FreeCounts,RPA_I,RPA_SC,RPA_DC,RAD51_Mon_I,RAD51_Mon_SC,RAD51_Mon_DC,RAD51_Dim_I,RAD51_Dim_SC,RAD51_Dim_DC] = LatticeSearch_Cluster(DNA,n_RAD51,n_A,n_D);
   [BoundCounts,RAD51_Mon_Bound,RAD51_Dim_Bound,RPA_A_Bound_MicroDiss,RPA_A_Open_MicroBind,RPA_D_Bound_MicroDiss,RPA_D_Open_MicroBind] = LatticeSearch_Bound(DNA,n_RAD51,n_A,n_D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MacroBind_Propensity = [k_on_RAD51*L_RAD51_M*FreeCounts(4),k_on_RAD51*L_RAD51_M*FreeCounts(5)*w_RAD51,k_on_RAD51*L_RAD51_M*FreeCounts(6)*(w_RAD51^2),...
                            k_on_RAD51*L_RAD51_D*FreeCounts(7),k_on_RAD51*L_RAD51_D*FreeCounts(8)*w_RAD51,k_on_RAD51*L_RAD51_D*FreeCounts(9)*(w_RAD51^2),...
                            k_on_RPA_A*L_RPA*FreeCounts(1),k_on_RPA_A*L_RPA*FreeCounts(2)*w_RPA,k_on_RPA_A*L_RPA*FreeCounts(3)];    %propensity functions of macroscopic binding
    MicroRPA_Propensity = [k_off_RPA_A*BoundCounts(3),k_off_RPA_D*BoundCounts(5),k_on_RPA_A*BoundCounts(4),k_on_RPA_D*BoundCounts(6)];  %propensity functions of RPA microscopic binding/unbinding
    Unbind_Propensity = [k_off_RAD51*BoundCounts(1),k_off_RAD51*BoundCounts(2)];    %propensity functions of RAD51 unbinding
    Full_Propensity = [MacroBind_Propensity,MicroRPA_Propensity,Unbind_Propensity]; %Full Propensity function (each one labeled below)
%         1: RAD51 Monomer Isolated Binding     7: RPA Macro Isolated Binding       %14: RAD51 Monomer Unbinding
%         2: RAD51 Monomer SC Binding           8: RPA Macro SC Binding             %15: RAD51 Dimer Unbinding
%         3: RAD51 Monomer DC Binding           9: RPA Macro DC Binding
%         4: RAD51 Dimer Isolated Binding       10: RPA-A Micro Unbinding
%         5: RAD51 Dimer SC Binding             11: RPA-D Micro Unbinding
%         6: RAD51 Dimer DC Binding             12: RPA-A Micro Binding
%                                               13: RPA-D Micro Binding
    
    a_0 = sum(Full_Propensity); %sum of all propensity functions
    Randoms = [rand,rand];  %random numbers for Monte Carlo step
    dt(Event_Count) = (1/a_0)*log(1/Randoms(1)); %time until next reaction occurs
    
    RPA = sort([RPA_I,RPA_SC,RPA_DC]);  %sorted list of all available RPA binding locations
    for S = 1:numel(Full_Propensity)
        SuM(S) = sum(Full_Propensity(1:S));
    end
    j = find(SuM >= (Randoms(2)*a_0),1,'first');  %reaction number, j, follows the same order as the propensity functions listed above
    ReactHist(Event_Count) = j;
    %Direct Method of Calculating which reaction is next to occur
    if j == 1   %RAD51 Monomer Isolated Binding
        Location = RAD51_Mon_I(randi(numel(RAD51_Mon_I)));
        DNA(2,Location:Location+(n_RAD51-1)) = RAD51;
    elseif j == 2   %RAD51 Monomer SC Binding
        Location = RAD51_Mon_SC(randi(numel(RAD51_Mon_SC)));
        DNA(2,Location:Location+(n_RAD51-1)) = RAD51;
    elseif j == 3   %RAD51 Monomer DC Binding
        Location = RAD51_Mon_DC(randi(numel(RAD51_Mon_DC)));
        DNA(2,Location:Location+(n_RAD51-1)) = RAD51;
    elseif j == 4   %RAD51 Dimer Isolated Binding
        Location = RAD51_Dim_I(randi(numel(RAD51_Dim_I)));
        DNA(2,Location:Location+(2*n_RAD51-1)) = RAD51;
    elseif j == 5   %RAD51 Dimer SC Binding
        Location = RAD51_Dim_SC(randi(numel(RAD51_Dim_SC)));
        DNA(2,Location:Location+(2*n_RAD51-1)) = RAD51;
    elseif j == 6   %RAD51 Dimer DC Binding
        Location = RAD51_Dim_DC(randi(numel(RAD51_Dim_DC)));
        DNA(2,Location:Location+(2*n_RAD51-1)) = RAD51;
    elseif j == 7   %RPA Macro Isolated Binding
        LocationA = RPA_I(randi(numel(RPA_I)));
        LocationD = LocationA + n_A;
        DNA(2,LocationA:LocationA+(n_A-1)) = RPA_A;
        DNA(2,LocationD:LocationD+(n_D-1)) = RPA_D;
    elseif j == 8   %RPA Macro SC Binding
        LocationA = RPA_SC(randi(numel(RPA_SC)));
        LocationD = LocationA + n_A;
        DNA(2,LocationA:LocationA+(n_A-1)) = RPA_A;
        DNA(2,LocationD:LocationD+(n_D-1)) = RPA_D;
    elseif j == 9   %RPA Macro DC Binding
        LocationA = RPA_DC(randi(numel(RPA_DC)));
        LocationD = LocationA + n_A;
        DNA(2,LocationA:LocationA+(n_A-1)) = RPA_A;
        DNA(2,LocationD:LocationD+(n_D-1)) = RPA_D;
    elseif j == 10  %RPA-A Micro Unbinding
        Location = RPA_A_Bound_MicroDiss(randi(numel(RPA_A_Bound_MicroDiss)));
        DNA(2,Location:Location+(n_A-1)) = 0;    DNA(1,Location:Location+(n_A-1)) = RPA_A;
        if DNA(1,Location+n_A) == RPA_D
            DNA(1,Location:Location+(n_RPA-1)) = 0; %full unbinding of RPA
        end
    elseif j == 11  %RPA-D Micro Unbinding
        Location = RPA_A_Bound_MicroDiss(randi(numel(RPA_A_Bound_MicroDiss)));
        DNA(2,Location:Location+(n_D-1)) = 0;    DNA(1,Location:Location+(n_D-1)) = RPA_A;
        if DNA(1,Location-1) == RPA_A
            DNA(1,Location-n_A:Location-n_A+(n_RPA-1)) = 0; %full unbinding of RPA
        end
    elseif j == 12  %RPA-A Micro Binding
        Location = RPA_A_Open_MicroBind(randi(numel(RPA_A_Open_MicroBind)));
        DNA(2,Location:Location+(n_A-1)) = RPA_A;   DNA(1,Location:Location+(n_A-1)) = 0;
    elseif j == 13  %RPA-D Micro Binding
        Location = RPA_D_Open_MicroBind(randi(numel(RPA_D_Open_MicroBind)));
        DNA(2,Location:Location+(n_D-1)) = RPA_D;   DNA(1,Location:Location+(n_D-1)) = 0;
    elseif j == 14  %RAD51 Monomer Unbindng
        Location = RAD51_Mon_Bound(randi(numel(RAD51_Mon_Bound)));
        DNA(2,Location:Location+(n_RAD51-1)) = 0;
    elseif j == 15  %RAD51 Dimer Unbinding
        Location = RAD51_Dim_Bound(randi(numel(RAD51_Dim_Bound)));
        DNA(2,Location:Location+(2*n_RAD51-1)) = 0;
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