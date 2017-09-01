clear all
close all


Rs=50;
RL=10;
L=0.1e-9;
C=10e-6;
w0=10e5;
Ton=1/w0;
Toff=1/w0;
dt = 0.005e-6;
Vsrc=1;

Tsim=8e-6;

%steady-state solution coeffs for driven mode
temp = 1/(j*w0*L)/(1/Rs+1/(j*w0*L)+j*w0*C);
A=abs(temp*Vsrc/Rs);
phi = angle(temp);

%params for homogeneous solution in charging mode
tau = Rs*C;
wr = 1/sqrt(L*C);
wd=0.5*sqrt(4*wr^2-(1/tau)^2);

%params for homogeneous solution in discharging mode (p is for "prime")
taup = RL*C;
wr = 1/sqrt(L*C);
wdp=0.5*sqrt(4*wr^2-(1/taup)^2);

t=0:dt:Tsim;
Tcycle=0;
IL=zeros(size(t));
V=zeros(size(t));

P=IL(1)-A*sin(w0*t(1)+phi);
Q=(1/wd)*(P/(2*tau)+V(1)/L-w0*cos(w0*t(1)+phi));


for ii=2:length(t)
    %update the local time        
    Tcycle=Tcycle+dt;
    if Tcycle<Ton %charging
        if Tcycle<dt
            %compute coefficients
            P=IL(ii-1)-A*sin(w0*t(ii)+phi);
            Q=(1/wd)*(P/(2*tau)+V(ii-1)/L-w0*A*cos(w0*t(ii)+phi));
        end
        %compute current and voltage in charging mode
        IL(ii)=exp(-Tcycle/(2*tau)).*(P*cos(wd*Tcycle)+Q*sin(wd*Tcycle))+A*sin(w0*t(ii)+phi);
        V(ii) = L*exp(-Tcycle/(2*tau)).*(P*(-wd*sin(wd*Tcycle)-1/(2*tau)*cos(wd*Tcycle))+...
                Q*(wd*cos(wd*Tcycle)-1/(2*tau)*sin(wd*Tcycle)))+...
                L*w0*A*cos(w0*t(ii)+phi);
    elseif Tcycle<=Ton+Toff %discharging
        if Tcycle<Ton+dt
            %compute coefficients
            d1=IL(ii-1);
            d2=(V(ii-1)/L+IL(ii-1)/(2*taup))/wdp;
        end
        %compute current and voltage in discharge mode
        IL(ii)=exp(-Tcycle/(2*taup)).*(d1*cos(wd*Tcycle)+d2*sin(wd*Tcycle));
        V(ii)=L*exp(-Tcycle/(2*taup)).*(d1*(-wdp*sin(wdp*Tcycle)-1/(2*taup)*cos(wdp*Tcycle))+...
              d2*(-wdp*cos(wdp*Tcycle)-1/(2*taup)*sin(wdp*Tcycle)));
    else
        %end of off period - reset local time to restart period (set to -dt
        %since the first thing we do is increment it)
        Tcycle=-dt;
    end
end

figure; 
plot(t,IL)
ylabel('Inductor current (A)')
xlabel('time (s)')
figure;
plot(t,V)
ylabel('LC Voltage (V)')
xlabel('time (s)')

fstep=1/Tsim;

temp = abs(fft(V));
Vf = temp(1:floor(length(temp)/2));


temp = abs(fft(IL));
If = temp(1:floor(length(temp)/2));


figure; 
plot((0:(length(Vf)-1))*fstep,Vf)
title('FFT of LC Voltage')
figure; 
plot((0:(length(If)-1))*fstep,If)
title('FFT of Inductor Current')