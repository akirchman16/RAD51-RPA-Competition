clear all;
close all;

% This code will model the competition between RAD51 and RPA in two stages.
% The first stage will include binding RPA to the lattice at random
% locations to act as an "initial condition" for the lattice. The second
% stage will actually model the competition between the hinging part of RPA
% and the RAD51 molecules. These two will be able to easily bind and unbind
% from the lattice (assuming there's a free location) but the entire RPA
% molecule will not unbind from the lattice. Both stages will include
% cooperativity in the model. No cooperativity is included for the d part
% of RPA to rebind since it has no choice in where to bind. Dummy zeros
% are included in the lattice to account for cooperation.


N = 8660;   %length of DNA lattice
DNA = zeros(1,N+2);
minIterations = 1000;    %minimum number of events that will occur in competition stage

% % % 1st STAGE % % % - Runs until maximum coverage of RPA
RPA_n = 20; %length of RPA molecule
RPA_BindProb = 0.8; %binding probability of RPA molecules

xRPA_Bound = 0;

RPAa_n = ceil(RPA_n*0.5); %size of A piece of RPA
RPAd_n = RPA_n-RPAa_n;  %size of D piece of RPA

Repetitions = 0;
RPA_Saturation = length(find(DNA(2:N+1)~= 0))/N;
RPA_Equilibrium = 0;

while ~RPA_Equilibrium
    RPA_BindSpot(Repetitions+1) = randi([2,(N+1)-(RPA_n-1)]);    %random location to attempt binding
    if DNA(RPA_BindSpot(Repetitions+1):RPA_BindSpot(Repetitions+1)+(RPA_n-1))==0 & rand <= RPA_BindProb  %checks for free location and binding probability
        DNA(RPA_BindSpot(Repetitions+1):RPA_BindSpot(Repetitions+1)+(RPAa_n-1)) = 1; %binds RPA to lattice (A piece)
        DNA(RPA_BindSpot(Repetitions+1)+(RPAa_n):RPA_BindSpot(Repetitions+1)+(RPAa_n)+(RPAd_n-1)) = 3;    %D piece of RPA
        xRPA_Bound = xRPA_Bound+1;
        RPA_BindHist(Repetitions+1) = RPA_BindSpot(Repetitions+1);
    end
    BoundRPA(Repetitions+1) = xRPA_Bound;  %counts how many RPA molecules are currently bound
    RPA_Saturation(Repetitions+1) = length(find(DNA(2:N+1)~= 0))/N;
    
    if Repetitions > 100
        RPA_SaturationStates = RPA_Saturation(Repetitions-99:Repetitions+1);
        RPA_SaturationChange = abs(diff(RPA_SaturationStates(1:10:end)));
        if mean(RPA_SaturationChange) <= RPA_n/N
            RPA_Equilibrium = 1;
        else
            RPA_Equilibrium = 0;
        end
    end
    
    Repetitions = Repetitions+1;
end
Final_RPA_Saturation = length(find(DNA(2:N+1)~=0))/N;  %final RPA fractional coverage for the system

x = 1:Repetitions;

% figure();     %plots RPA Saturation with each iteration
% scatter(x,RPA_Saturation,5,'r','filled');
% xlabel('Iterations');
% ylabel('RPA Fractional Coverage');
% title('RPA Coverage vs. Iterations');

% % % 2nd STAGE % % %
% BoundAtSpot arrays
RPAd_BoundAtSpot = zeros(1,N);
RAD51_BoundAtSpot = zeros(1,N);

% RAD51 Properties and Arrays
RAD51_n = 3;    %length of protein
RAD51_k_on = 1; %kinetic rate constant for binding
RAD51_k_off = 1;    %kinetic rate constant for unbinding
RAD51_L = 1;    %concentration of molecules
w = 1;  %cooperativity parameter

RAD51_K = RAD51_k_on/RAD51_k_off;   %equilibrium constant

% RPA Properties and Arrays
RPAd_n = RPA_n-RPAa_n;    %size of D piece of RPA (part that hinges open)
RPAd_k_on = 1;  %kinetic rate constant for binding
RPAd_k_off = 1; %kinetic rate constant for unbinding

RPAd_K = RPAd_k_on/RPAd_k_off;  %equilibrium constant

% memory allocation arrays
a_RPAdUnbind = zeros(1,minIterations);
a_RAD51Bind_I = zeros(1,minIterations);
a_RAD51Unbind = zeros(1,minIterations);
a_RPAdBind = zeros(1,minIterations);
a_0 = zeros(1,minIterations);
xRPAd_Bound = zeros(1,minIterations); %number of bound molecules for each protein
xRAD51_Bound = zeros(1,minIterations);
xRPAd_Free = zeros(1,minIterations); %number of free locations for each protein
xRAD51_Free = zeros(1,minIterations);
t = zeros(1,minIterations);
dt = zeros(1,minIterations);
j = zeros(1,minIterations);
RPAd_BindSpot = zeros(1,minIterations);
RPAd_UnbindSpot = zeros(1,minIterations);
RAD51_BindSpot = zeros(1,minIterations);
RAD51_UnbindSpot = zeros(1,minIterations);
RAD51_FracCover = zeros(1,minIterations);
RAD51_BoundAtSpot = zeros(1,N+2);
Total_FracCover = zeros(1,minIterations);

% initial conditions for the lattice (already saturated with RPA - D parts all bound)
xRPAd_Bound(1) = xRPA_Bound;    %number of bound D-pieces of RPA
xRAD51_Bound(1) = 0;    %number of bound RAD51 proteins
xRPAd_Free(1) = 0;  %how many locations where RPAd can rebind
xRAD51_Free(1) = 0; %how many locations where RAD51 can bind
t(1) = 0;
RPAd_BindCounter = 0;
RPAd_UnbindCounter = 0;
RAD51_BindCounter = 0;
RAD51_BindCounter_I = 0;
RAD51_BindCounter_SC = 0;
RAD51_BindCounter_DC = 0;
RAD51_UnbindCounter = 0;
xRPAd_Free(1) = 0;
xRAD51_Free(1) = 0;
RAD51_FracCover(1) = 0;

% Preparing for the actual reactions
Events = 1;
Equilibrium = 0;

while ~Equilibrium  %loops until the system reaches equilibrium
%     counts how many of each RAD51 free spots there are
    RAD51_Free = 0;
    Free_Counter = 0;
    Isolated = 0;
    I_Counter = 0;
    SinglyContiguous = 0;
    SC_Counter = 0;
    DoublyContiguous = 0;
    DC_Counter = 0;
    for b = 2:(N+1)-(RAD51_n-1)
        if DNA(b:b+(RAD51_n-1)) == 6  %first checks if location is free for RAD51 binding
            if (DNA(b-1) == 1 && DNA(b+RAD51_n) == 6) || (DNA(b-1) == 6 && (DNA(b+RAD51_n) == 6 || DNA(b+RAD51_n) == 1 || DNA(b+RAD51_n) == 0))
                I_Counter = I_Counter+1;
                Isolated(Events,I_Counter) = b;
            elseif (DNA(b-1) == 10 && (DNA(b+RAD51_n) == 6 || DNA(b+RAD51_n) == 1 || DNA(b+RAD51_n) == 0)) || ((DNA(b-1) == 6 || DNA(b-1) == 1) && DNA(b+RAD51_n) == 10)
                SC_Counter = SC_Counter+1;
                SinglyContiguous(Events,SC_Counter) = b;
            elseif (DNA(b-1) == 10 && DNA(b+RAD51_n) == 10)
                DC_Counter = DC_Counter+1;
                DoublyContiguous(Events,DC_Counter) = b;
            end
            Free_Counter = Free_Counter+1;
            RAD51_Free(Free_Counter) = b;
        end
    end
    xRAD51_Free_I(Events) = I_Counter;
    xRAD51_Free_SC(Events) = SC_Counter;
    xRAD51_Free_DC(Events) = DC_Counter;
    
%     propensity functions/probabilities of reactions
    a_RPAdUnbind(Events) = RPAd_k_off*xRPAd_Bound(Events);  %RPA hinges open
    a_RAD51Bind_I(Events) = RAD51_k_on*RAD51_L*xRAD51_Free_I(Events);   %RAD51 binds - isolated
    a_RAD51Bind_SC(Events) = RAD51_k_on*RAD51_L*xRAD51_Free_SC(Events)*w; %RAD51 binds - singly contiguous
    a_RAD51Bind_DC(Events) = RAD51_k_on*RAD51_L*xRAD51_Free_DC(Events)*(w^2);   %RAD51 binds - doubly contiguous
    a_RAD51Unbind(Events) = RAD51_k_off*xRAD51_Bound(Events);   %RAD51 unbinds
    a_RPAdBind(Events) = RPAd_k_on*xRPAd_Free(Events);  %RPA hinges closed
    a_0(Events) = a_RPAdUnbind(Events)+a_RAD51Bind_I(Events)+a_RAD51Bind_SC(Events)+a_RAD51Bind_DC(Events)+a_RAD51Unbind(Events)+a_RPAdBind(Events);
    
%     First Reaction Method to choose what reaction occurs
    R_1 = rand;
    R_2 = rand;
    R_3 = rand;
    R_4 = rand;
    R_5 = rand;
    R_6 = rand;
    tau = [(1/a_RPAdUnbind(Events))*log(1/R_1),(1/a_RAD51Bind_I(Events))*log(1/R_2),(1/a_RAD51Bind_SC(Events))*log(1/R_3),(1/a_RAD51Bind_DC(Events))*log(1/R_4),(1/a_RAD51Unbind(Events))*log(1/R_5),(1/a_RPAdBind(Events))*log(1/R_6)];
    dt(Events) = min(tau);
    j(Events) = find(tau == min(tau));
    
%     Stores where proteins are bound
    RPAd_BoundAtSpot = zeros(1,N);  %clears RPAd_BoundAtSpot
    RPAd_OnSpots = find(diff([1 DNA 1]) == 2);   %list of where RPAd is bound
    RPAd_BoundAtSpot(RPAd_OnSpots) = 1;  %all locations where RPAd is bound
    RAD51_LeftEnds = find(diff([1 DNA 1]) == 5 | diff([1 DNA 1] == -4));   %left ends of all available RAD51 binding locations
    RPAd_LeftEnds = find(diff([1 DNA 1]) == 5);   %left ends of all possible binding locations for RPAd
    
%     Follow through with whichever reaction was chosen to occur
    if j(Events) == 1    %D-part of RPA unbinds from the lattice
        pos_RPAd_Unbind = randi(length(RPAd_OnSpots));
        RPAd_UnbindSpot(Events) = RPAd_OnSpots(pos_RPAd_Unbind);    %random location for RPA to unbind is determined
        DNA(RPAd_UnbindSpot(Events):RPAd_UnbindSpot(Events)+(RPAd_n-1)) = 6;   %unbinds D-piece of RPA
        
%         Update populations for each type of location/spot
        xRPAd_Bound(Events+1) = xRPAd_Bound(Events)-1;
        xRAD51_Bound(Events+1) = xRAD51_Bound(Events);
        xRPAd_Free(Events+1) = sum(max(find(diff([1 DNA 1]) == -6 | diff([1 DNA 1]) == -5 | diff([1 DNA 1]) == 4)-find(diff([1 DNA 1]) == 5 | diff([1 DNA 1]) == -4)-RPAd_n+1,0));   %number of free spots which RPAd can rebind

%         Bug checking
        RPAd_UnbindCounter = RPAd_UnbindCounter+1;     
    elseif j(Events) == 2    %RAD-51 protein binds to the lattice - isolated
        pos_RAD51_Bind = randi(length(Isolated(Events,:)));
        RAD51_BindSpot(Events) = Isolated(Events,pos_RAD51_Bind);
        
%         Binding happens
        DNA(RAD51_BindSpot(Events):RAD51_BindSpot(Events)+(RAD51_n-1)) = 10;
        RAD51_BoundAtSpot(RAD51_BindSpot(Events)) = 1;  %shows RAD51 is now bound at location

%         Update populations for each type of location/spot
        xRPAd_Bound(Events+1) = xRPAd_Bound(Events);
        xRAD51_Bound(Events+1) = xRAD51_Bound(Events)+1;
        xRPAd_Free(Events+1) = sum(max(find(diff([1 DNA 1]) == -6 | diff([1 DNA 1]) == -5 | diff([1 DNA 1]) == 4)-find(diff([1 DNA 1]) == 5 | diff([1 DNA 1]) == -4)-RPAd_n+1,0));   %number of free spots which RPAd can rebind

%         Bug Checking
        RAD51_BindCounter = RAD51_BindCounter+1;
        RAD51_BindCounter_I = RAD51_BindCounter_I+1;
    elseif j(Events) == 3    %RAD-51 protein binds to the lattice - singly contiguous
        pos_RAD51_Bind = randi(length(SinglyContiguous(Events,:)));
        RAD51_BindSpot(Events) = SinglyContiguous(Events,pos_RAD51_Bind);
        
%         Binding happens
        DNA(RAD51_BindSpot(Events):RAD51_BindSpot(Events)+(RAD51_n-1)) = 10;
        RAD51_BoundAtSpot(RAD51_BindSpot(Events)) = 1;  %shows RAD51 is now bound at location

%         Update populations for each type of location/spot
        xRPAd_Bound(Events+1) = xRPAd_Bound(Events);
        xRAD51_Bound(Events+1) = xRAD51_Bound(Events)+1;
        xRPAd_Free(Events+1) = sum(max(find(diff([1 DNA 1]) == -6 | diff([1 DNA 1]) == -5 | diff([1 DNA 1]) == 4)-find(diff([1 DNA 1]) == 5 | diff([1 DNA 1]) == -4)-RPAd_n+1,0));   %number of free spots which RPAd can rebind

%         Bug Checking
        RAD51_BindCounter = RAD51_BindCounter+1;
        RAD51_BindCounter_SC = RAD51_BindCounter_SC+1;
    elseif j(Events) == 4    %RAD-51 protein binds to the lattice - doubly contiguous
        pos_RAD51_Bind = randi(length(DoublyContiguous(Events,:)));
        RAD51_BindSpot(Events) = DoublyContiguous(Events,pos_RAD51_Bind);
        
%         Binding happens
        DNA(RAD51_BindSpot(Events):RAD51_BindSpot(Events)+(RAD51_n-1)) = 10;
        RAD51_BoundAtSpot(RAD51_BindSpot(Events)) = 1;  %shows RAD51 is now bound at location

%         Update populations for each type of location/spot
        xRPAd_Bound(Events+1) = xRPAd_Bound(Events);
        xRAD51_Bound(Events+1) = xRAD51_Bound(Events)+1;
        xRPAd_Free(Events+1) = sum(max(find(diff([1 DNA 1]) == -6 | diff([1 DNA 1]) == -5 | diff([1 DNA 1]) == 4)-find(diff([1 DNA 1]) == 5 | diff([1 DNA 1]) == -4)-RPAd_n+1,0));   %number of free spots which RPAd can rebind

%         Bug Checking
        RAD51_BindCounter = RAD51_BindCounter+1;
        RAD51_BindCounter_DC = RAD51_BindCounter_DC+1;
    elseif j(Events) == 5    %RAD-51 protein unbinds from the lattice
        RAD51_OnSpots = find(RAD51_BoundAtSpot == 1);
        pos_RAD51_Unbind = randi(length(RAD51_OnSpots));
        RAD51_UnbindSpot(Events) = RAD51_OnSpots(pos_RAD51_Unbind);   %random location for RAD51 to unbind
        
%         Unbinding happens
        DNA(RAD51_UnbindSpot(Events):RAD51_UnbindSpot(Events)+(RAD51_n-1)) = 6; %unbinds RAD51 from lattice
        RAD51_BoundAtSpot(RAD51_UnbindSpot(Events)) = 0;    %shows RAD51 is no longer bound at spot
        
%         Update populations for each type of location/spot
        xRPAd_Bound(Events+1) = xRPAd_Bound(Events);
        xRAD51_Bound(Events+1) = xRAD51_Bound(Events)-1;
        xRPAd_Free(Events+1) = sum(max(find(diff([1 DNA 1]) == -6 | diff([1 DNA 1]) == -5 | diff([1 DNA 1]) == 4)-find(diff([1 DNA 1]) == 5 | diff([1 DNA 1]) == -4)-RPAd_n+1,0));   %number of free spots which RPAd can rebind     
        
%         Bug Checking
        RAD51_UnbindCounter = RAD51_UnbindCounter+1;
    elseif j(Events) == 6    %D-part of RPA re-binds to the lattice
        RPAd_FreeCounter = 0;
        RPAd_FreeSpots = 0;
        for d = 1:length(RPAd_LeftEnds)
            if DNA(RPAd_LeftEnds(d):RPAd_LeftEnds(d)+(RPAd_n-1)) == 6    %checks if each possible position is free
                RPAd_FreeCounter = RPAd_FreeCounter+1;
                RPAd_FreeSpots(RPAd_FreeCounter) = RPAd_LeftEnds(d);   %all possilble locations for RAD51 to bind
            end
        end
        pos_RPAd_Bind = randi(length(RPAd_FreeSpots));
        RPAd_BindSpot(Events) = RPAd_FreeSpots(pos_RPAd_Bind);

%         Binding happens
        DNA(RPAd_BindSpot(Events):RPAd_BindSpot(Events)+(RPAd_n-1)) = 3;    %rebinds RPAd to lattice

%         Update populations for each type of location/spot
        xRPAd_Bound(Events+1) = xRPAd_Bound(Events)+1;
        xRAD51_Bound(Events+1) = xRAD51_Bound(Events);
        xRPAd_Free(Events+1) = sum(max(find(diff([1 DNA 1]) == -6 | diff([1 DNA 1]) == -5 | diff([1 DNA 1]) == 4)-find(diff([1 DNA 1]) == 5 | diff([1 DNA 1]) == -4)-RPAd_n+1,0));   %number of free spots which RPAd can rebind

%          Bug Checking
        RPAd_BindCounter = RPAd_BindCounter+1;
    end
    
    t(Events+1) = t(Events)+dt(Events); %advance time
    RAD51_FracCover(Events+1) = length(find(DNA == 10))/N;  %calculates RAD51 saturation
    Total_FracCover(Events+1) = length(find(DNA ~= 0))/N;   %calculates total saturation of the DNA
    
    if Events >= minIterations   %loop to check if system is in equilibrium or not
        PreviousStates = RAD51_FracCover(ceil(0.25*Events):Events-1);   %states of the system over that las 100 events
        EqTest = PreviousStates(1:10:end);    %state of the system each of the last 10 events for the last 100 events
        RAD51_FracCover_Change = max(PreviousStates)-min(PreviousStates); %changes of the state for the last 10 events the last 100 events
        if (mean(abs(diff(EqTest))) <= 2*RAD51_n/N) || Events >= 10000
            Equilibrium = 1;
        else
            Equilibrium = 0;
        end
    end
    
    Events = Events+1;
end

RAD51_EqFracCover = mean(RAD51_FracCover(Events-100:Events));    %Equilibrium Fractional Coverage of RAD51

% Calculates how many times each type of event occured in the simulation
RPAd_UnbindingEvents = length(find(j == 1));
RAD51_BindingEvents = length(find(j == 2))+length(find(j == 3))+length(find(j == 4));
RAD51_BindingEvents_I = length(find(j == 2));
RAD51_BindingEvents_SC = length(find(j == 3));
RAD51_BindingEvents_DC = length(find(j == 4));
RAD51_UnbindingEvents = length(find(j == 5));
RPAd_BindingEvents = length(find(j == 6));
EventTypes = [RPAd_UnbindingEvents RAD51_BindingEvents_I RAD51_BindingEvents_SC RAD51_BindingEvents_DC RAD51_UnbindingEvents RPAd_BindingEvents];

figure('Name','Bounded Exponential Fit');
scatter(t,RAD51_FracCover,3,'r','filled');
hold on;
xlabel('Time, t (s)');
xlim([0 max(t)]);
legend({'Data','Fit'},'Location','east');
ylabel('RAD-51 Fractional Coverage');
title('RAD-51 Fractional Coverage');