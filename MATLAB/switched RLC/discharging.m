function [V, IL] = discharging( V0, I0L, RL, L, C, Toff)
%Calculates voltage and current after "discharging" phase of switched parallel
%RLC resonator, given circuit parameters and ending voltage and inductor
%current of previous discharge cycle. Careful: w0 is the driving frequency 
%not the resonant frequency. t is the start time of the "on" period.
%See Kim & Wang, "Theory of Switched RF Resonators", IEEE Trans. Circuits
%and Systems, v. 53, no. 12, 2006.

% %params for homogeneous solution in charging mode
% tau = Rs*C;
% wr = 1/sqrt(L*C);
% gam1 = -1/(2*tau)+0.5*sqrt((1/tau)^2-4*wr^2);
% gam2 = -1/(2*tau)-0.5*sqrt((1/tau)^2-4*wr^2);

%params for homogeneous solution in discharging mode (p is for "prime")
taup = RL*C;
wr = 1/sqrt(L*C);
wdp=0.5*sqrt(4*wr^2-(1/taup)^2);

%calculate coefficients for solution
d1=I0L;
d2=(V0/L+I0L/(2*taup))/wdp;

V=L*d1*exp(-Toff/(2*taup))*(-wdp*sin(wdp*Toff)-1/(2*taup)*cos(wdp*Toff))+...
    L*d2*exp(-Toff/(2*taup))*(wdp*cos(wdp*Toff)-1/(2*taup)*sin(wdp*Toff));

IL=exp(-Toff/(2*taup))*(d1*cos(wdp*Toff)+d2*sin(wdp*Toff));

end

