function [ V, IL ] = charging( V0, I0L, Rs, L, C, Ton, Vsrc, w0, t)
%Calculates voltage and current after "charging" phase of switched parallel
%RLC resonator, given circuit parameters and ending voltage and inductor
%current of previous discharge cycle. Careful: w0 is the driving frequency 
%not the resonant frequency. t is the start time of the "on" period.
%See Kim & Wang, "Theory of Switched RF Resonators", IEEE Trans. Circuits
%and Systems, v. 53, no. 12, 2006.

%amplitude and phase of steady-state part of inductor current (particular solution)
temp = 1/(j*w0*L)/(1/Rs+1/(j*w0*L)+j*w0*C);
A=abs(temp*Vsrc/Rs);
phi = angle(temp);

%params for homogeneous solution in charging mode
tau = Rs*C;
wr = 1/sqrt(L*C);
wd=0.5*sqrt(4*wr^2-(1/tau)^2);

%calculate coefficients for solution
P=I0L-A*sin(phi);
Q=(1/wd)*(P/(2*tau)+V0/L-w0*cos(w0*t+phi));

%calculate IL and V at end of "on" period
IL = P*exp(-Ton/(2*tau))*cos(wd*Ton)+Q*exp(-Ton/(2*tau))*sin(wd*Ton)+A*sin(w0*Ton+w0*t+phi);

V = P*exp(-Ton/(2*tau))*(-wd*sin(wd*Ton)-1/(2*tau)*cos(wd*Ton))+...
    Q*exp(-Ton/(2*tau))*(wd*cos(wd*Ton)-1/(2*tau)*sin(wd*Ton))+...
    A*w0*cos(w0*t+w0*Ton+phi);

end

