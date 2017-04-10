%%%% BASIC SIR GILLESPIE ALGORITHM
%%%% This script is a basic stochastic SIR model for measles
%%%% It assumes the events occur exponentially and generates times
%%%% accordingly.
%%%% It caluclates the probability of each event occuring and
%%%% selects the event. 
%%%% It outputs a graph of the mean number of susceptibles, infecteds and
%%%% recovereds over time after N simulations. 
%% Parameters
bet = 800; %infection rate
gamm = 365/7; %recovery rate
M = 10000; %Max time
P = 5000; %Population
N = 100; %Number of Simualtions

%% Initialise Matrices
span = 0:0.001:0.05; %Grid to map individiual paths onto
S_sum = zeros(length(span), 1); %Vector for Rolling Sum
S_interp = zeros(length(span), N); %Vector for interpolated values
I_sum = zeros(length(span), 1);
I_interp = zeros(length(span), N);
R_sum = zeros(length(span), 1);
R_interp = zeros(length(span), N);
a = zeros(1, 2); %empty propensity vector

%% Gillespie Algorithm
for k = 1:N
    rr = rand(M, 2); %generate random numbers
    t = zeros(M, 1); %empty time vector
    SIR = zeros(M, 3); %Empty Classes for M months for each of the 3 classes
    SIR(1, :)= [0.95*P 0.05*P 0]; %intial sizes for each class 0.1% of the pop are infected
    i = 1;
    while t(i) < max(span)
        %Compute propensity functions
        a(1) = SIR(i, 1)*SIR(i, 2)*bet/P;
        a(2) = SIR(i, 2)*gamm;
        a0 = a(1)+a(2);  %compute total propensity
        t(i+1) = t(i) + (1/a0)*log(1/rr(i, 1)); %time until next change
        if k == M
            return
        end
        %Update Classes
        if rr(i, 2) < sum(a(1))/a0
            %If random number is less than recovery proability
            %an infection occurs...
            SIR(i+1, 1) = SIR(i, 1) - 1; %An S becomes I
            SIR(i+1, 2) = SIR(i, 2) + 1;
            SIR(i+1, 3) = SIR(i, 3);
        else
            %Else a recovery occurs
            SIR(i+1, 1) = SIR(i, 1); %An I becomes R
            SIR(i+1, 2) = SIR(i, 2) - 1;
            SIR(i+1, 3) = SIR(i, 3) + 1;
        end
    i = i + 1;
    end
    %Interpolate onto a grid of uniform time steps
    %and calculate moving sum
    S_interp(:, k) = interp1(t(1:i), SIR(1:i, 1), span, 'previous');
    S_sum = S_sum + S_interp(:, k);
    I_interp(:, k) = interp1(t(1:i), SIR(1:i, 2), span, 'previous');
    I_sum = I_sum + I_interp(:, k);
    R_interp(:, k) = interp1(t(1:i), SIR(1:i, 3), span, 'previous');
    R_sum = R_sum + R_interp(:, k);
end
S_mean = S_sum/N; %Calculate Mean
I_mean = I_sum/N;
R_mean = R_sum/N;


%% Plots
%Plot Mean
plot(span, S_mean, span, I_mean, span, R_mean)
xlabel('time')
ylabel('Number of Individuals')
legend('Susceptibles', 'Infecteds', 'Recvoreds');