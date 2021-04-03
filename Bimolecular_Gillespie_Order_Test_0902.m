clear all;
close all;

% This code is set to test whether the subscripts of reactions matter, if
% the order in which the propensity functions are tested, matters.

% Inputs that are the same with either test

Omega = 1;     %volume of container in which reaction occurs
k_f = 0.1;    %reaction rate constant for the forward reaction (A+B-->C)
k_r = 0.1;    %reaction rate constant for the reverse reaction (C-->A+B)

c_f = k_f/Omega;    %calculation of c for forward reaction
c_r = k_r;          %calculation of c for reverse reaction

% Runs the model by testing the forward reaction propensity first and then
% the sum of the two. The second part will test the reverse reaction
% propensity function before the sum of the two.

xA_1(1) = 1000;   %initial population of A
xB_1(1) = 500;   %initial population of B
xC_1(1) = 0;     %initial population of C

t_1(1) = 0;           %initializes time at t=0

ODExA_1(1) = xA_1(1);   %intial populations for ODE approximation
ODExB_1(1) = xB_1(1);
ODExC_1(1) = xC_1(1);

Equilibrium_1 = 0;
Loops_1 = 0;

ForwardCount_1 = 0;
ReverseCount_1 = 0;

while ~Equilibrium_1
    a_f_1(Loops_1+1) = c_f*xA_1(Loops_1+1)*xB_1(Loops_1+1);   %propensity functions
    a_r_1(Loops_1+1) = c_r*xC_1(Loops_1+1);
    a_0_1(Loops_1+1) = a_f_1(Loops_1+1)+a_r_1(Loops_1+1);     %sum of propensity functions
    
    r_1_1(Loops_1+1) = rand; %random numbers to be used to calculate tau and j
    r_2_1(Loops_1+1) = rand;
    
    tau_1(Loops_1+1) = (1/a_0_1(Loops_1+1))*log(1/r_1_1(Loops_1+1));  %Monte Carlo method to find tau
    if a_f_1(Loops_1+1) > r_2_1(Loops_1+1)*a_0_1(Loops_1+1)  %determines which reaction occurs
        j_1(Loops_1+1) = 1;           %forward reaction
    else
        j_1(Loops_1+1) = 2;           %reverse reaction
    end
    
    if j_1(Loops_1+1) == 1               %updates populations for a forward reaction
        xA_1((Loops_1+1)+1) = xA_1(Loops_1+1)-1;
        xB_1((Loops_1+1)+1) = xB_1(Loops_1+1)-1;
        xC_1((Loops_1+1)+1) = xC_1(Loops_1+1)+1;
        ForwardCount_1 = ForwardCount_1+1;
    elseif j_1(Loops_1+1) == 2           %updates populations for a reverse reaction
        xA_1((Loops_1+1)+1) = xA_1(Loops_1+1)+1;
        xB_1((Loops_1+1)+1) = xB_1(Loops_1+1)+1;
        xC_1((Loops_1+1)+1) = xC_1(Loops_1+1)-1;
        ReverseCount_1 = ReverseCount_1+1;
    end
    t_1((Loops_1+1)+1) = t_1(Loops_1+1)+tau_1(Loops_1+1);   %advances time by time step, tau
    
    state_1(Loops_1+1,1) = xA_1(Loops_1+1);   %creates history of state vectors after each
    state_1(Loops_1+1,2) = xB_1(Loops_1+1);   %reaction
    state_1(Loops_1+1,3) = xC_1(Loops_1+1);
    
    dx_1 = [(-c_f*tau_1(Loops_1+1)),(c_r*tau_1(Loops_1+1)); (-c_f*tau_1(Loops_1+1)),(c_r*tau_1(Loops_1+1)); (c_f*tau_1(Loops_1+1)),(-c_r*tau_1(Loops_1+1))]*[ODExA_1(Loops_1+1)*ODExB_1(Loops_1+1);ODExC_1(Loops_1+1)];
    ODEold_1 = [ODExA_1(Loops_1+1); ODExB_1(Loops_1+1); ODExC_1(Loops_1+1)];    %Euler method to solve ODEs
    ODEnew_1 = dx_1 + ODEold_1;
    ODExA_1((Loops_1+1)+1) = ODEnew_1(1,1);
    ODExB_1((Loops_1+1)+1) = ODEnew_1(2,1);
    ODExC_1((Loops_1+1)+1) = ODEnew_1(3,1);
    
    ODEstate_1(Loops_1+1,1) = ODExA_1(Loops_1+1);
    ODEstate_1(Loops_1+1,2) = ODExB_1(Loops_1+1);
    ODEstate_1(Loops_1+1,3) = ODExC_1(Loops_1+1);
    
    if Loops_1 > 100
        xC_EqTest_1 = [xC_1(Loops_1-100) xC_1(Loops_1)]; %state of system 50 loops ago and now
        xC_Change_1 = diff(xC_EqTest_1);    %difference between the two states
        if abs((xC_Change_1)) <= 1    %if system hasn't changed by more than 1 molecule in 50 iterations it's at equilibrium
            Equilibrium_1 = 1;
        else
            Equilibrium_1 = 0;
        end
    end
    
    Loops_1 = Loops_1+1;
end

eqxA_1 = round(mean(xA_1(round(0.75*Loops_1):Loops_1)));  %equilibrium pop.
eqxB_1 = round(mean(xB_1(round(0.75*Loops_1):Loops_1)));
eqxC_1 = round(mean(xC_1(round(0.75*Loops_1):Loops_1)));
eqODExA_1 = round(mean(ODExA_1(round(0.75*Loops_1):Loops_1)));   %ODE equilibrium
eqODExB_1 = round(mean(ODExB_1(round(0.75*Loops_1):Loops_1)));
eqODExC_1 = round(mean(ODExC_1(round(0.75*Loops_1):Loops_1)));

% Tests a run of testing reverse propensity function first in order to
% compare to testing the forward one first.

xA_2(1) = xA_1(1);   %initial population of A
xB_2(1) = xB_1(1);   %initial population of B
xC_2(1) = xC_1(1);     %initial population of C

t_2(1) = 0;           %initializes time at t=0

ODExA_2(1) = xA_2(1);   %intial populations for ODE approximation
ODExB_2(1) = xB_2(1);
ODExC_2(1) = xC_2(1);

Equilibrium_2 = 0;
Loops_2 = 0;

ForwardCount_2 = 0;
ReverseCount_2 = 0;

while ~Equilibrium_2
    a_f_2(Loops_2+1) = c_f*xA_2(Loops_2+1)*xB_2(Loops_2+1);   %propensity functions
    a_r_2(Loops_2+1) = c_r*xC_2(Loops_2+1);
    a_0_2(Loops_2+1) = a_f_2(Loops_2+1)+a_r_2(Loops_2+1);     %sum of propensity functions
    
    r_1_2(Loops_2+1) = rand; %random numbers to be used to calculate tau and j
    r_2_2(Loops_2+1) = rand;
    
    tau_2(Loops_2+1) = (1/a_0_2(Loops_2+1))*log(1/r_1_2(Loops_2+1));  %Monte Carlo method to find tau
    if a_r_2(Loops_2+1) > r_2_2(Loops_2+1)*a_0_2(Loops_2+1)  %determines which reaction occurs
        j_2(Loops_2+1) = 2;           %reverse reaction
    else
        j_2(Loops_2+1) = 1;           %forward reaction
    end
    
    if j_2(Loops_2+1) == 1               %updates populations for a forward reaction
        xA_2((Loops_2+1)+1) = xA_2(Loops_2+1)-1;
        xB_2((Loops_2+1)+1) = xB_2(Loops_2+1)-1;
        xC_2((Loops_2+1)+1) = xC_2(Loops_2+1)+1;
        ForwardCount_2 = ForwardCount_2+1;
    elseif j_2(Loops_2+1) == 2           %updates populations for a reverse reaction
        xA_2((Loops_2+1)+1) = xA_2(Loops_2+1)+1;
        xB_2((Loops_2+1)+1) = xB_2(Loops_2+1)+1;
        xC_2((Loops_2+1)+1) = xC_2(Loops_2+1)-1;
        ReverseCount_2 = ReverseCount_2+1;
    end
    t_2((Loops_2+1)+1) = t_2(Loops_2+1)+tau_2(Loops_2+1);   %advances time by time step, tau
    
    state_2(Loops_2+1,1) = xA_2(Loops_2+1);   %creates history of state vectors after each
    state_2(Loops_2+1,2) = xB_2(Loops_2+1);   %reaction
    state_2(Loops_2+1,3) = xC_2(Loops_2+1);
    
    dx_2 = [(-c_f*tau_2(Loops_2+1)),(c_r*tau_2(Loops_2+1)); (-c_f*tau_2(Loops_2+1)),(c_r*tau_2(Loops_2+1)); (c_f*tau_2(Loops_2+1)),(-c_r*tau_2(Loops_2+1))]*[ODExA_2(Loops_2+1)*ODExB_2(Loops_2+1);ODExC_2(Loops_2+1)];
    ODEold_2 = [ODExA_2(Loops_2+1); ODExB_2(Loops_2+1); ODExC_2(Loops_2+1)];    %Euler method to solve ODEs
    ODEnew_2 = dx_2 + ODEold_2;
    ODExA_2((Loops_2+1)+1) = ODEnew_2(1,1);
    ODExB_2((Loops_2+1)+1) = ODEnew_2(2,1);
    ODExC_2((Loops_2+1)+1) = ODEnew_2(3,1);
    
    ODEstate_2(Loops_2+1,1) = ODExA_2(Loops_2+1);
    ODEstate_2(Loops_2+1,2) = ODExB_2(Loops_2+1);
    ODEstate_2(Loops_2+1,3) = ODExC_2(Loops_2+1);
    
    if Loops_2 > 100
        xC_EqTest_2 = [xC_2(Loops_2-100) xC_2(Loops_2)]; %state of system 50 loops ago and now
        xC_Change_2 = diff(xC_EqTest_2);    %difference between the two states
        if abs((xC_Change_2)) <= 1    %if system hasn't changed by more than 1 molecule in 50 iterations it's at equilibrium
            Equilibrium_2 = 1;
        else
            Equilibrium_2 = 0;
        end
    end
    
    Loops_2 = Loops_2+1;
end

eqxA_2 = round(mean(xA_2(round(0.75*Loops_2):Loops_2)));  %equilibrium pop.
eqxB_2 = round(mean(xB_2(round(0.75*Loops_2):Loops_2)));
eqxC_2 = round(mean(xC_2(round(0.75*Loops_2):Loops_2)));
eqODExA_2 = round(mean(ODExA_2(round(0.75*Loops_2):Loops_2)));   %ODE equilibrium
eqODExB_2 = round(mean(ODExB_2(round(0.75*Loops_2):Loops_2)));
eqODExC_2 = round(mean(ODExC_2(round(0.75*Loops_2):Loops_2)));

figure();                          %plot of molecular populations
scatter(t_1,xA_1,5,'red','filled');
hold on;
scatter(t_1,xB_1,5,'green','filled');
scatter(t_1,xC_1,5,'blue','filled');
plot(t_1,ODExA_1,'red');
plot(t_1,ODExB_1,'green');
plot(t_1,ODExC_1,'blue');
scatter(t_2,xA_2,5,'yellow','filled');
scatter(t_2,xB_2,5,'cyan','filled');
scatter(t_2,xC_2,5,'magenta','filled');
plot(t_2,ODExA_2,'yellow');
plot(t_2,ODExB_2,'cyan');
plot(t_2,ODExC_2,'magenta');
xlabel('Time');
ylabel('Populations');
ylim([0 max(max(state_1))]);
title('A+B<-->C');
legend('xA','xB','xC');