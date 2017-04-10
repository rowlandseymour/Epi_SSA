function [ dxdt ] = SIR_ODE(t, x, b, bet, gamm, mu)
%The correspnding ODE function for Measles ODE
% dS/dt = b -muS - betaSI, dI/dt = betaSI - muI - gammaI, dR/dt = gammaI
% - muR
dxdt(1, 1) = b - mu*x(1, 1) - bet*x(1, 1)*x(2, 1);
dxdt(2, 1) = bet*x(1, 1)*x(2, 1) - mu*x(2, 1) - gamm*x(2, 1);
dxdt(3, 1) = gamm*x(2, 1) - mu*x(3, 1);
end