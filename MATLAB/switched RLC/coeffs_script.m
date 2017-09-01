clear all
close all

Rs=50;
RL=10;
L=1e-9;
C=10e-12;
Ton=1e-4;
Toff=1e-4;
Vsrc=1;
w0=1e6;
N=100;


IL = zeros(1,N);
V=IL;

IL(1)=0;
V(1)=0;
[V(2), IL(2)] = charging(0,0,Rs, L, C, Ton, Vsrc, w0, 0);
[V(3), IL(3)] = discharging(V(2), IL(2), RL, L, C, Ton);

for ii=2:floor(N/2)
    [V(2*(ii)), IL(2*(ii))] = charging(V(2*(ii)-1), IL(2*(ii)-1),Rs, L, C, Ton, Vsrc, w0, (ii-1)*(Ton+Toff));
    [V(2*(ii)+1), IL(2*(ii)+1)] = discharging(V(2*(ii)), IL(2*(ii)), RL, L, C, (ii-1)*(Ton+Toff)+Ton);
end

clock = repmat([1 1 0 0],1,floor(N/2));

figure;
plot(1:length(V), V, (1:0.5:(N+0.5)), clock*max(abs(V)))
legend({'V';'clock'})

figure;
plot(1:length(IL), IL, 1:0.5:(N+0.5), clock*max(abs(IL)))
legend({'I_L';'clock'})