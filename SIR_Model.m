%% INFECTED -> RECOVERY MODEL
%This script solves the differential equation dI/dt - -gamma * I(t)
%It comparse this to the mean of N simulations from an SSA. 

%% Set Parameters and Matrices
N = 100; %Number of Simulations
M = 25; %Number of Time Steps
n0 = 25; %Initial Number of Infecteds
gamm = 1; %Rate of Recoverey
r = rand(M, N); %Generate Random Numbers
I = repmat([n0:-1:0].', 1, N); %Set number of infecteds at each time point
t = zeros(M, N); %initalise time vector

%% SSA
%Calculate the time at which each infected recovers
for k = 1:N
    for i = 1:25 %Infecteds cannot go negative
        t(i + 1, k) = t(i, k) + 1/(I(i, k)*gamm)*log(1/r(i, k));
    end
end

%% Interpolation
%Interpolate this onto uniform grid so we can calculate mean
span = 0:0.01:6; %Uniform Grid of time points
I_interp = zeros(length(span), N);
for i = 1:N
    I_interp(:,i) = interp1(t(:, i), I(:, i), span, 'previous'); %interpolate onto gird
end
I_interp(isnan(I_interp)) = 0; %Set outside points to 0
I_mean = sum(I_interp, 2)/N; %mean number of infecteds

%% Soving the ODE
[t_ode, x] = ode45(@(t, x) -gamm*x, span, n0); %Call solver

%% Plots
%Plot 10 paths from SSA simulations
%Compare ODE to mean of SSA simulations
plot(span, I_interp(:, 1:10)) %Plot 10 Paths
xlabel('t')
ylabel('I(t)')
figure
hold on
plot(span, I_mean) %Plot Mean
plot(t_ode, x) %Plot ODE Solution
hold off
legend('Mean Number of Infecteds', 'ODE Solution')
xlabel('t')
ylabel('Number of Infecteds')