clear all
close all

Rs=50;
RL=10;
L=1e-6;
C=22e-9;
alpha=0.5; %duty cycle (less than 1)
vs=1;
is=vs/Rs;
w0=1/sqrt(L*C);
wc=2*pi*100e3;
tau = C/(1/Rs+1/RL);
wd = sqrt(w0^2-(1/(2*tau))^2);

T=1/(2e6);
dt=T/100;

%specify simulation time
N=200;

mon=floor(alpha*T/dt);
moff=floor((1-alpha)*T/dt);
m=mon+moff;

t=0:dt:(N*m*dt);
ton = dt*(0:(mon-1));
toff = dt*(0:(moff-1));

%particular (steady-state) solution coefficient

temp = 1/(1-wc^2*L*C+j*wc*L/Rs+j*wc*L/RL);
A=abs(temp);
phi = angle(temp);

iL = zeros(size(t));
v = zeros(size(t));

B1p=zeros(1,N);
B2p=zeros(1,N);
B1m=zeros(1,N);
B2m=zeros(1,N);

for n=1:N
    %charging phase (+ state)
    %set constants
    if n==1 %assume initially iL=0 and v=0
        B1p(n) = -A*is*sin(phi);
        B2p(n) = (1/wd)*(A*is*sin(phi)/(2*tau)-wc*A*is*cos(phi));
    else
        B1p(n) = iL((n-1)*m)-A*is*sin(wc*(n-1)*T+phi);
        B2p(n) = (1/wd)*(v((n-1)*m)/L-B1p(n)/(2*tau)-...
                wc*A*is*cos(wc*(n-1)*T+phi));
    end
    %fill iL and v for this "on" chunk of time
    iL((1:mon)+(n-1)*m)=exp(-ton/(2*tau)).*(B1p(n)*cos(wd*ton)+B2p(n)*sin(wd*ton))...
        +A*is*sin(wc*t((1:mon)+(n-1)*m)+phi);
    v((1:mon)+(n-1)*m)=-L*exp(-ton/(2*tau))/(2*tau).*(B1p(n)*cos(wd*ton)+B2p(n)*sin(wd*ton))+...
        L*exp(-ton/(2*tau)).*(-wd*B1p(n)*sin(wd*ton)+wd*B2p(n)*cos(wd*ton))+...
        L*wc*A*is*cos(wc*t((1:mon)+(n-1)*m)+phi);
    
    %switch flips, changing the sign of Vs
    
    %discharging phase (- state)
    %set constants
    B1m(n) = iL((n-1)*m+mon)+A*is*sin(wc*(n-1+alpha)*T+phi);
    B2m(n) = (1/wd)*(v((n-1)*m+mon)/L-B1m(n)/(2*tau)+...
        wc*A*is*cos(wc*(n-1+alpha)*T+phi));
    
    iL((n-1)*m+mon+(1:moff))=exp(-toff/(2*tau)).*(B1m(n)*cos(wd*toff)+B2m(n)*sin(wd*toff))...
        -A*is*sin(wc*t((n-1)*m+mon+(1:moff))+phi);
    v((n-1)*m+mon+(1:moff))=-L*exp(-toff/(2*tau))/(2*tau).*(B1m(n)*cos(wd*toff)+B2m(n)*sin(wd*toff))+...
        L*exp(-toff/(2*tau)).*(-wd*B1m(n)*sin(wd*toff)+wd*B2m(n)*cos(wd*ton))-...
        L*wc*A*is*cos(wc*t((n-1)*m+mon+(1:moff))+phi);
end
figure
plot(t, iL)
xlabel('Time (s)')
ylabel('Inductor Current (A)')

figure
plot(t, v)
xlabel('Time (s)')
ylabel('Resonator Voltage (V)')

figure
plot((1:N)*T, B1p, (1:N)*T, B2p)
xlabel('Time (s)')
ylabel('Solution constants (+ state)')

figure
plot((1:N)*T, B1m, (1:N)*T, B2m)
xlabel('Time (s)')
ylabel('Solution constants (- state)')