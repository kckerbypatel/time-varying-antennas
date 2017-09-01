function [ V, IL, clock ] = coeffs( V0, IL0, Rs, RL, L, C, Ton, Toff, Vsrc, w0, N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
IL = zeros(1,N);
V=IL;

IL(1)=0;
V(1)=0;
[V(2), IL(2)] = charging(V0,IL0,Rs, L, C, Ton, Vsrc, w0, 0);
[V(3), IL(3)] = discharging(V(2), IL(2), RL, L, C, Ton);

for ii=2:floor(N/2)
    [V(2*(ii)), IL(2*(ii))] = charging(V(2*(ii)-1), IL(2*(ii)-1),Rs, L, C, Ton, Vsrc, w0, (ii-1)*(Ton+Toff));
    [V(2*(ii)+1), IL(2*(ii)+1)] = discharging(V(2*(ii)), IL(2*(ii)), RL, L, C, (ii-1)*(Ton+Toff)+Ton);
end

clock = repmat([1 1 0 0],1,floor(N/2));

end

