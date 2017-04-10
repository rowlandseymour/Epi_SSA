function [ dxdt ] = SQUIRREL_ODE(t, x, a, bet, sig, alph, K, b, Boxes)
%SQUIRRELPOX ODE
% Solves the squirrelpox matrix ODE with mixing and disperseal
%Create Mixing Matrix
alph_main = (1-alph)*ones(Boxes, 1);
alph_1 = alph/2*ones(Boxes, 1);
alph_mat = spdiags([alph_1,alph_main,alph_1], -1:1, Boxes,Boxes);
alph_mat(1, 2) = alph;
alph_mat(Boxes, Boxes-1) = alph;
%Use for barrier example
% alph_mat(9, 8) = 0.9*alph/2;
% alph_mat(9, 10) = 0.1*alph/2;
% alph_mat(10, 9) = 0.1*alph/2;
% alph_mat(10, 11) = 0.1*alph/2;
% alph_mat(11, 10) = 0.1*alph/2;
% alph_mat(11, 12) = 0.9*alph/2;
%Vectorise Parameters
a = a*ones(Boxes, 1);
bet = bet*ones(Boxes, 1);
sig = sig*ones(Boxes, 1);
b = b*ones(Boxes,1);
q = (a-b)./K;
%Change Notation to S and I
S = x(1:Boxes);
I = x(Boxes+1 :2*Boxes);
N = S+I;
%Calculate Dispersal Coefficients
D = zeros(Boxes, 1);
for i = 1:Boxes
    if i == 1 || i == Boxes %Special Case for 1st and Last Box
        D(i) = 0.5*b(i)*exp((N(i)-K(i))/(0.5*K(i)));
    else
        D(i) = b(i)*exp((N(i)-K(i))/(0.5*K(i)));
    end
end
d_mat = spdiags([D/2,-D,D/2], -1:1, Boxes,Boxes);
d_mat(1, 1) = -0.5*D(1); %Set End Cases
d_mat(Boxes, Boxes) = -0.5*D(end);
%d mat(10, 10) = 0;
%Set up Matrix ODE
BSI = bet.*S.*(alph_mat*I);
dxdt = [(a - (q.*N)).*N - b.*S- BSI + d_mat*S; BSI - (b + sig).*I+d_mat*I];
end