%%% SQUIRREL POX MODEL
%%% This scripts simulates an outbreak of squirrel pox on an island of red
%%% squirrels. 
%%% It partions the islands into a 1D row of boxes, and allows the quirrels
%%% to mix with the neighbouring boxes. 
%%% It simualates N simualtions fo the Gillepsie algorithm with birth,
%%% natural death, pox-induced death, and transmission. 
%%% The boxes iniaitlly contain the same numeber of squirrels and have the
%%% same carrying capacity.
%%% The squirrels are also allowed to disperse, permenantly moving to a
%%% non-neighbouring box. 
%%% This scirpt also solves the ODE system for this model and plots the
%%% results as travelling waves. 
%% Parameters
Boxes = 50; %Number of Boxes
P = 19*ones(Boxes, 1); %Initial Population
K = 20*ones(Boxes, 1); %Carrying Capacity
aa = 1.5; %Birth
death = 0.9; %Death
d = 0.9; %Dispersal
bet = 5; %Transmission
sig = 26; %Pox Death
q = (aa-death)./K; %Crowding
I_initial = zeros(Boxes, 1); %Only Box 1 infected
I_initial(1, 1) = 3; %3 Grey Squirrels Introduced
p = 0.3; %mixing rate

%% Initialise Matrices and Variable
span = 0:0.001:10; %Grid to map individiual paths onto
N = 10; %number of simulations
M = 100000; %Max time
S = zeros(M, Boxes); %Susceptible Matrix
I = zeros(M, Boxes); %Infected Matrix
S_interp = zeros(length(span), Boxes); %Uniform Susceptible Matrix
I_interp = zeros(length(span), Boxes); %Uniform Infected Matrix
S_mean = zeros(length(span), Boxes); %Mean Number of Suscpetibles
I_mean = zeros(length(span), Boxes); %Mean Number of Infecteds
a = zeros(7, Boxes); %Propensity vector
count = zeros(7, N); %Count how many times each event occurs
s_count = zeros(4);
i_count = zeros(4);
S(1, :) = P; %Set intial Conditions
I(1, :) = I_initial;


%% Gillespie Algorithm
for k = 1:N
    rr = rand(M, 2); %generate random numbers
    t = zeros(M, 1); %empty time vector
    i = 1;
    while t(i) < max(span)
        for b = 1:Boxes
            H = S(i, b) + I(i, b); %Total Population
            if b == 1 || b == Boxes %Special End Cases
                D = 0.5*d*exp(-(K(b)-H)/(0.5*K(b)));
            else
                D = d*exp(-(K(b)-H)/(0.5*K(b)));
            end
        %Calcualte Propensity for each box
        a(1, b) = ((aa-q(b)*H)*H);  %Growth
        a(2, b) = death*S(i, b);    %S Death
        a(3, b) = death*I(i, b);    %I Death
        a(4, b) = sig*I(i, b);      %I Pox Death
        a(5, b) = D*S(i, b);        %S Dispersal
        a(6, b) = D*I(i, b);        %I Dispersal
        %Infection
        if b == 1
            a(7, b) = bet*S(i, b)*(p*I(i, 2)+I(i, 1));
        elseif b == Boxes
            a(7, b) = bet*S(i, b)*(p*I(i, Boxes-1)+I(i, b));
        else
            a(7, b) = bet*S(i, b)*(p/2*(I(i, b-1)+ I(i, b+1))+I(i, b));
        end
        end
        a0 = sum(a, 1); %Total Propensity
        ca0 = cumsum(a0); %Cummualtive Propensity
        ii = find(ca0/ca0(Boxes) >= rr(i, 2)); % Find the Box of the event
        box = ii(1);
        r_local = -1*(rr(i, 2) - ca0(box)/ca0(Boxes)); %Transfrom r
        t(i+1) = t(i) - (1/ca0(Boxes))*log(rr(i, 1)); %time at next change
        %Gillespie Algorithm
            S(i+1, :) = S(i, :); %Carry over old populations
            I(i+1, :) = I(i, :);
        if r_local < a(1, box)/ca0(Boxes)
            S(i+1, box) = S(i, box) + 1; %Birth
        count(1) = count(1) + 1;
        elseif r_local < sum(a(1:2, box))/ca0(Boxes)
            S(i+1, box) = S(i, box) - 1; %S Death
            count(2) = count(2) + 1;
        elseif r_local < sum(a(1:3, box))/ca0(Boxes)
            I(i+1, box) = I(i, box) - 1; %Death of I
            count(3) = count(3) + 1;
        elseif r_local < sum(a(1:4, box))/ca0(Boxes)
            I(i+1, box) = I(i, box) - 1; %pox induced Death
            count(4) = count(4) + 1;
        elseif r_local < sum(a(1:5, box))/ca0(Boxes)
            r_local2 = ca0(Boxes)/a(5, box)*(sum(a(1:5, box))/ca0(Boxes) - r_local);
            count(5) = count(5) + 1;
        if box == 1
            s_count(1) = s_count(1)+1;
            S(i+1, 2) = S(i, 2)+1; %Dispersal of S
            S(i+1, 1) = S(i, 1)-1;
        elseif box == Boxes
            s_count(2) = s_count(2)+1;
            S(i+1, box-1) = S(i, box-1)+1; %Dispersal of S
            S(i+1, box) = S(i, box)-1;
        elseif r_local2 < 0.5
            s_count(3) = s_count(3)+1;
            S(i+1, box-1) = S(i, box-1)+1; %Dispersal of S
            S(i+1, box) = S(i, box)-1;
        else
            s_count(4) = s_count(4)+1;
            S(i+1, box+1) = S(i, box+1)+1; %Dispersal of S
            S(i+1, box) = S(i, box)-1;
        end
        elseif r_local < sum(a(1:6, box))/ca0(Boxes)
            r_local3 = ca0(Boxes)/a(5, box)*(sum(a(1:5, box))/ca0(Boxes) - r_local);
            count(6) = count(6) + 1;
        if box == 1
            i_count(1) = i_count(1)+1;
            I(i+1, box+1) = I(i, box+1)+1; %Dispersal of I
            I(i+1, box) = I(i, box)-1;
        elseif box == Boxes
            i_count(2) = i_count(2)+1;
            I(i+1, box-1) = I(i, box-1)+1; %Dispersal of I
            I(i+1, box) = I(i, box)-1;
        elseif r_local3 < 0.5
            i_count(3) = i_count(3)+1;
            I(i+1, box-1) = I(i, box-1)+1; %Dispersal of I
            I(i+1, box) = I(i, box)-1;
        else
            i_count(4) = i_count(4)+1;
            I(i+1, box+1) = I(i, box+1)+1; %Dispersal of I
            I(i+1, box) = I(i, box)-1;
        end
        else
            count(7) = count(7) + 1;
            S(i+1, box) = S(i, box) - 1; %Infection
            I(i+1, box) = I(i, box) + 1;
        end
            i = i+1;
        end
        %Interpolate the popualtions onto a unifrom time grid
    for b = 1:Boxes
        S_interp(:, b) = interp1(t(1:i), S(1:i, b), span, 'previous');
        I_interp(:, b) = interp1(t(1:i), I(1:i, b), span, 'previous');
    end
    S_mean = S_mean + S_interp/N; %calculate mean
    I_mean = I_mean + I_interp/N;
end
%% ODE Solution
x0 = [P; I_initial]; %Initial Conditions
s_mat = spdiags(ones(Boxes, 3), -1:1, Boxes, Boxes); %Jacoabian Pattern
s_mat = kron([1 1 ; 1 1], s_mat);
options = odeset('JPattern', s_mat, 'reltol', 1.e-10, 'abstol', 1.e-40);
[t_ode, x] = ode15s(@SQUIRREL_ODE, 0:0.001:20, x0, options, aa, bet, sig, 0.3, K, death, Boxes);
%% Plots
plot(span, S_mean) %Plot Mean
figure
plot(span, I_mean)
xlabel('Time (Years)')
ylabel('Mean Number of Infecteds')
%Calculate Travelling Wave
S_tw= zeros(Boxes, 1);
for i = 1:Boxes
S_tw(i) = min(span(S_mean(:, i)<17).'); %First time S(t)<17
end
figure
plot(S_tw)
xlabel('Cell')
ylabel('min t S(t)<17')