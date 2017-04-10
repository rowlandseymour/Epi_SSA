%%% SQUIRREL POX MODEL
%%% This scripts simulates an outbreak of squirrel pox on an island of red
%%% squirrels. 
%%% It partions the islands into a 1D row of boxes, and allows the quirrels
%%% to mix with the neighbouring boxes. 
%%% It simualates N simualtions fo the Gillepsie algorithm with birth,
%%% natural death, pox-induced death, and transmission. 
%%% The boxes iniaitlly contain the same numeber of squirrels and have the
%%% same carrying capacity.
%%% This model introuces small-world connections, where certain boxes
%%% non-neighbouring are connected, allowing long-rnage dispersal between
%%% them. 
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
sw_connections = [28 3 24 42 26 21 43 8 35 12];
%% Initialise Matrices and Variable
span = 0:0.00001:25; %Grid to map individiual paths onto
N = 1; %number of simulations
M = 100000; %Number of Time Steps
S = zeros(M, Boxes); %Susceptible Matrix
I = zeros(M, Boxes); %Infected Matrix
S_interp = zeros(length(span), Boxes); %Uniform Susceptible Matrix
I_interp = zeros(length(span), Boxes); %Uniform Infected Matrix
S_mean = zeros(length(span), Boxes); %Mean Number of Suscpetibles
I_mean = zeros(length(span), Boxes); %Mean Number of Infecteds
a = zeros(7, Boxes); %Propensity vector
count = zeros(7, N); %Count how many times each event occurs
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
        %Calcualte Propensity
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
    ii = find(ca0/ca0(Boxes) >= rr(i, 2)); % Find the Box
    box = ii(1);
    r_local = -1*(rr(i, 2) - ca0(box)/ca0(Boxes)); %Transform r
    t(i+1) = t(i) - (1/ca0(Boxes))*log(rr(i, 1)); %time at next change
    %Gillespie Algorithm
        S(i+1, :) = S(i, :); %Carry over old populations
        I(i+1, :) = I(i, :);
    if r_local < a(1, box)/ca0(Boxes)
        S(i+1, box) = S(i, box) + 1; %Birth
        count(1, k) = count(1, k) + 1;
    elseif r_local < sum(a(1:2, box))/ca0(Boxes)
        S(i+1, box) = S(i, box) - 1; %S Death
        count(2, k) = count(2, k) + 1;
    elseif r_local < sum(a(1:3, box))/ca0(Boxes)
    	I(i+1, box) = I(i, box) - 1; %Death of I
        count(3, k) = count(3, k) + 1;
    elseif r_local < sum(a(1:4, box))/ca0(Boxes)
        I(i+1, box) = I(i, box) - 1; %pox induced Death
        count(4, k) = count(4, k) + 1;
    elseif r_local < sum(a(1:5, box))/ca0(Boxes)%Dispersal of S
        count(5, k) = count(5, k) + 1;
    %Transform Random Number
    r_local2 = ca0(Boxes)/a(5, box)*(sum(a(1:5, box))/ca0(Boxes) - r_local);
    %If the box has no long-range connections, then choose
    %which of the neighbours it disperses to.
    if ismember(box, sw_connections)==0
        if box == 1
            S(i+1, box+1) = S(i, box+1)+1;
            S(i+1, box) = S(i, box)-1;
        elseif box == Boxes
            S(i+1, box-1) = S(i, box-1)+1;
            S(i+1, box) = S(i, box)-1;
        elseif r_local2 < 0.5
            S(i+1, box-1) = S(i, box-1)+1;
            S(i+1, box) = S(i, box)-1;
        else
            S(i+1, box+1) = S(i, box+1)+1;
            S(i+1, box) = S(i, box)-1;
        end
    elseif ismember(box, sw_connections(1:2:end))==1
    connection = find(sw_connections == box);
    %If the box is in the first half of the pairs of
    %numbers, then chose whether to disperse to a neighbour
    %or long-range partner
    if box == 1
    %First and last boxes either move to only neighbour
    %or long-dispersal
        if r_local2 < 0.5
            S(i+1, box+1) = S(i, box+1)+1;
            S(i+1, box) = S(i, box)-1;
        else
            S(i+1, sw_connections(connection+1))= S(i, sw_connections(connection+1))+1;
            S(i+1, box) = S(i, box)-1;
        end
    elseif box == Boxes
        if r_local2 < 0.5
            S(i+1, box-1) = S(i, box-1)+1;
            S(i+1, box) = S(i, box)-1;
        else
            S(i+1, sw_connections(connection+1))= S(i, sw_connections(connection+1))+1;
            S(i+1, box) = S(i, box)-1;
        end
    %Other cells move to either of the neighbours or
    %the long-range partner
    elseif r_local2 < 1/3
        S(i+1, box-1) = S(i, box-1)+1; %Dispersal of S
        S(i+1, box) = S(i, box)-1;
    elseif r_local2 < 2/3
        S(i+1, box+1) = S(i, box+1)+1; %Dispersal of S
        S(i+1, box) = S(i, box)-1;
    else
        S(i+1, sw_connections(connection+1))= S(i, sw_connections(connection+1))+1;
        S(i+1, box) = S(i, box)-1;
    end
    else
    %If the box is in the second half of the pairs of
    %numbers, then chose whether to disperse to a neighbour
    %or long-range partner
    connection = find(sw_connections == box);
    if box == 1
        if r_local2 < 0.5
            S(i+1, box+1) = S(i, box+1)+1; %Dispersal of S
            S(i+1, box) = S(i, box)-1;
        else
            S(i+1, sw_connections(connection-1))= S(i, sw_connections(connection-1))+1;
            S(i+1, box) = S(i, box)-1;
        end
    elseif box == Boxes
        if r_local2 < 0.5
            S(i+1, box-1) = S(i, box-1)+1; %Dispersal of S
            S(i+1, box) = S(i, box)-1;
        else
            S(i+1, sw_connections(connection-1))= S(i, sw_connections(connection-1))+1;
            S(i+1, box) = S(i, box)-1;
        end
    elseif r_local2 < 1/3
        S(i+1, box-1) = S(i, box-1)+1; %Dispersal of S
        S(i+1, box) = S(i, box)-1;
    elseif r_local2 < 2/3
        S(i+1, box+1) = S(i, box+1)+1; %Dispersal of S
        S(i+1, box) = S(i, box)-1;
    else
        S(i+1, sw_connections(connection-1))= S(i, sw_connections(connection-1))+1;
        S(i+1, box) = S(i, box)-1;
    end
    end
    elseif r_local < sum(a(1:6, box))/ca0(Boxes) % I Dispersal
    %This event mirrors the S Dispersal Event
    r_local3 = ca0(Boxes)/a(6, box)*(sum(a(1:6, box))/ca0(Boxes) - r_local);
    count(6, k) = count(6, k) + 1;
    if ismember(box, sw_connections)==0
        if box == 1
            I(i+1, box+1) = I(i, box+1)+1;
            I(i+1, box) = I(i, box)-1;
        elseif box == Boxes
            I(i+1, box-1) = I(i, box-1)+1;
            I(i+1, box) = I(i, box)-1;
        elseif r_local3 < 0.5
            I(i+1, box-1) = I(i, box-1)+1;
            I(i+1, box) = I(i, box)-1;
        else
            I(i+1, box+1) = I(i, box+1)+1;
            I(i+1, box) = I(i, box)-1;
        end
    elseif ismember(box, sw_connections(1:2:end))==1
        connection = find(sw_connections == box);
    if box == 1
        if r_local3 < 0.5
            I(i+1, box+1) = I(i, box+1)+1;
            I(i+1, box) = I(i, box)-1;
        else
            I(i+1, sw_connections(connection+1))= I(i, sw_connections(connection+1))+1;
            I(i+1, box) = I(i, box)-1;
        end
    elseif box == Boxes
    if r_local3 < 0.5
        I(i+1, box-1) = I(i, box-1)+1;
        I(i+1, box) = I(i, box)-1;
    else
        I(i+1, sw_connections(connection+1))= I(i, sw_connections(connection+1))+1;
        I(i+1, box) = I(i, box)-1;
    end
    elseif r_local3 < 1/3
        I(i+1, box-1) = I(i, box-1)+1;
        I(i+1, box) = I(i, box)-1;
    elseif r_local3 < 2/3
        I(i+1, box+1) = I(i, box+1)+1;
        I(i+1, box) = I(i, box)-1;
    else
        I(i+1, sw_connections(connection+1))= I(i, sw_connections(connection+1))+1;
        I(i+1, box) = I(i, box)-1;
    end
    else
        connection = find(sw_connections == box);
    if box == 1
        if r_local3 < 0.5
            I(i+1, box+1) = I(i, box+1)+1;
            I(i+1, box) = I(i, box)-1;
        else
            I(i+1, sw_connections(connection-1))= I(i, sw_connections(connection-1))+1;
            I(i+1, box) = I(i, box)-1;
        end
    elseif box == Boxes
        if r_local3 < 0.5
            I(i+1, box-1) = I(i, box-1)+1;
            I(i+1, box) = I(i, box)-1;
        else
        I(i+1, sw_connections(connection-1))= I(i, sw_connections(connection-1))+1;
            I(i+1, box) = I(i, box)-1;
    end
    elseif r_local3 < 1/3
        I(i+1, box-1) = I(i, box-1)+1;
        I(i+1, box) = I(i, box)-1;
    elseif r_local3 < 2/3
        I(i+1, box+1) = I(i, box+1)+1;
        I(i+1, box) = I(i, box)-1;
    else
            I(i+1, sw_connections(connection-1)) = I(i, sw_connections(connection-1))+1;
            I(i+1, box) = I(i, box)-1;
        end
    end
    else
    count(7, k) = count(7, k) + 1;
        S(i+1, box) = S(i, box) - 1; %Infection
        I(i+1, box) = I(i, box) + 1;
    end
        i = i+1;
    end
    for b = 1:Boxes %Interpolate S and I onto a uniform time grid
        S_Interp(:, b) = interp1(t(1:i), S(1:i, b), span, 'previous');
        I_Interp(:, b) = interp1(t(1:i), I(1:i, b), span, 'previous');
    end
%Calculate Mean
S_mean = S_mean + S_Interp/N;
I_mean = I_mean + I_Interp/N;
end
%% Plots
plot(span, S_mean) %Plot Mean
figure
plot(span, I_mean)
xlabel('Time (Years)')
ylabel('Mean Number of Infecteds')
beep
%Travelling Wave Plots
S_tw = zeros(Boxes, 1);
for i = 1:Boxes
S_tw(i) = min(span(S_mean(:, i)<17).'); %Calculate First Time S(t)<17
end
S_tw_connections = zeros(Boxes, 1);
for i = 1:Boxes
if ismember(i, sw_connections) == 0
S_tw_connections(i) = 0;
else
S_tw_connections(i) = i;
end
end
S_tw_connections(S_tw_connections == 0)= NaN;