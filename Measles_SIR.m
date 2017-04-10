%%% MEASLES SIR MODEL
%%% This model uses the Gillespie algorithm to simulate a measles outbreak
%%% in a large city.
%%% It simulates N Gillespie Simulations for an SIR Model with infection
%%% rate beta, recovery rate gamma, and recovery rate mu. 
%%% It compares this to the equivalent ODE Model.
%%% It is used to show the ODE tends to a stationary point, bu the
%%% stochastic model oscilates. 

%% Initialise Parameters
clear all
close all
P =1e6; %Size of Population
bet = 782.3304; %infection rate
mu = 1/80; %death rate
b = mu*P; %birth rate
gamm = 365/7; %recovery rate


%% Initialise Matrices and Variables
span = 0:0.001:5; %Grid to map individiual paths onto
N = 50; %number of simulations
M = 5000000; %Number of Time Steps
S_sum = zeros(length(span), 1); %Empty vector for susceptibles sum
S_interp = zeros(length(span), N); %Empty vector for number of susceptibles
I_sum = zeros(length(span), 1);
I_interp = zeros(length(span), N);
R_sum = zeros(length(span), 1);
R_interp = zeros(length(span), N);
a = zeros(1, 6); %empty propensity vector


%% Gillespie Algorithm
for k = 1:N
rr = rand(M, 2); %generate random numbers
t = zeros(M, 1); %empty time vector
SIR = zeros(M, 3); %Empty Classes for M months for each of the 3 classes
SIR(1, :)= [0.065*P 80 P-0.065*P-80]; %initial sizes for each class 5
i = 1;
while t(i) < max(span)
    %Compute propensity functions
    a(1) = b;
    a(2) = SIR(i, 1)*SIR(i, 2)*bet/P;
    a(3) = SIR(i, 1)*mu;
    a(4) = SIR(i, 2)*gamm;
    a(5) = SIR(i, 2)*mu;
    a(6) = SIR(i, 3)*mu;
    a0 = sum(a);
    t(i+1) = t(i) + (1/a0)*log(1/rr(i, 1)); %time at next change
    if i == M
        fprintf('Time Steps Exceeded')
    end
    %Gillespie Algorithm
    if rr(i, 2) < a(1)/a0
        SIR(i+1, 1) = SIR(i, 1) + 1; %Birth
        SIR(i+1, 2) = SIR(i, 2);
        SIR(i+1, 3) = SIR(i, 3);
    elseif rr(i, 2) < sum(a(1:2))/a0
        SIR(i+1, 1) = SIR(i, 1) - 1; %Infection
        SIR(i+1, 2) = SIR(i, 2) + 1;
        SIR(i+1, 3) = SIR(i, 3);
    elseif rr(i, 2) < sum(a(1:3))/a0
        SIR(i+1, 1) = SIR(i, 1) - 1; %Death of S
        SIR(i+1, 2) = SIR(i, 2);
        SIR(i+1, 3) = SIR(i, 3);
    elseif rr(i, 2) < sum(a(1:4))/a0
        SIR(i+1, 1) = SIR(i, 1); %Recovery
        SIR(i+1, 2) = SIR(i, 2) - 1;
        SIR(i+1, 3) = SIR(i, 3) + 1;
    elseif rr(i, 2) < sum(a(1:5))/a0
        SIR(i+1, 1) = SIR(i, 1); %Death of I
        SIR(i+1, 2) = SIR(i, 2) - 1;
        SIR(i+1, 3) = SIR(i, 3);
    else
        SIR(i+1, 1) = SIR(i, 1); %Death of R
        SIR(i+1, 2) = SIR(i, 2);
        SIR(i+1, 3) = SIR(i, 3) - 1;
    end
        i = i+1;
    end
    S_interp(:, k) = interp1(t(1:i), SIR(1:i, 1), span, 'previous');%Interpolate
    %onto uniform grid
    S_sum = S_sum + S_interp(:, k); %calculate moving sum
    I_interp(:, k) = interp1(t(1:i), SIR(1:i, 2), span, 'previous');
    I_sum = I_sum + I_interp(:, k); %calculate moving sum
    R_interp(:, k) = interp1(t(1:i), SIR(1:i, 3), span, 'previous');
    R_sum = R_sum + R_interp(:, k); %calculate moving sum
end
S_mean = S_sum/N; %Calculate Mean
I_mean = I_sum/N; %Calculate Mean
R_mean = R_sum/N; %Calculate Mean


%% Solving the ODE Model
[t_ode, x] = ode23s(@SIR_ODE, span, SIR(1, :), [], b, bet/P, gamm, mu); %call solver


%% Plots
%Plot 10 Sample Paths from the SSA simulations
%Plot Mean of SSAs and ODE
plot(span, I_interp(:, 1:10));
xlabel('Time (years)')
ylabel('I(t)')
title('10 Simulations of Infected Individuals')
figure
plot(span, I_mean) %Plot Mean
hold on
plot(t_ode, x(:, 2))
hold off
xlabel('Time (years)')
ylabel('I(t)')
legend('Mean Number of Infecteds', 'ODE')