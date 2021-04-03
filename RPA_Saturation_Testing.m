clear all;
close all;

N = 1000;   %length of DNA lattice
DNA = zeros(1,N);
minIterations = 1000;    %minimum number of events that will occur in competition stage

% % % 1st STAGE % % % - Runs until maximum coverage of RPA
RPA_n = 30; %length of RPA molecule
RPA_BindProb = 0.5; %binding probability of RPA molecules

xRPA_Bound = 0;

RPAa_n = ceil(RPA_n/2); %size of A piece of RPA
RPAd_n = RPA_n-RPAa_n;  %size of D piece of RPA

Repetitions = 0;
RPA_Saturation = length(find(DNA~= 0))/N;
RPA_Equilibrium = 0;

while RPA_Equilibrium == 0
    RPA_BindSpot(Repetitions+1) = randi(N-(RPA_n-1));    %random location to attempt binding
    if DNA(RPA_BindSpot(Repetitions+1):RPA_BindSpot(Repetitions+1)+(RPA_n-1))==0 &  rand <= RPA_BindProb  %checks for free location and binding probability
        DNA(RPA_BindSpot(Repetitions+1):RPA_BindSpot(Repetitions+1)+(RPAa_n-1)) = 1; %binds RPA to lattice (A piece)
        DNA(RPA_BindSpot(Repetitions+1)+(RPAa_n):RPA_BindSpot(Repetitions+1)+(RPAa_n)+(RPAd_n-1)) = 3;    %D piece of RPA
        xRPA_Bound = xRPA_Bound+1;
        RPA_BindHist(Repetitions+1) = RPA_BindSpot(Repetitions+1);
    end
    BoundRPA(Repetitions+1) = xRPA_Bound;  %counts how many RPA molecules are currently bound
    RPA_Saturation(Repetitions+1) = length(find(DNA~= 0))/N;
    Repetitions = Repetitions+1;
    
    if Repetitions > 100    %equilibrium test
        if abs(RPA_Saturation(Repetitions-100)-RPA_Saturation(Repetitions)) <= (RPA_n/N)
            RPA_Equilibrium = 1;
        else
            RPA_Equilibrium = 0;
        end
    end
end
Final_RPA_Saturation = length(find(DNA~=0))/N;

x = 1:Repetitions;

figure();
scatter(x,RPA_Saturation,5,'r','filled');
xlabel('Iterations');
ylabel('RPA Fractional Coverage');
title('RPA Coverage vs. Iterations');