%% Uppgift 1 a)

% An = 0
% Bn = 4/(k*pi)
% Se uppgift1.txt
%% Uppgift 1 b)
clf
syms n x k

T=1;
w=2*pi/T;
M=200;
y=0;
f=0;
t=T*(0:M-1)/M;

for n=1:100
    An = 0
    Bn = 4*mod(n,2)/(n*pi) %mod(n,2) för att få bort jämna k
    y = y + An*cos(n*w*t) + Bn*sin(n*w*t);
end
plot(t,y)
hold on

%% Uppgift 2 a)
clf
clc

%num = (s.^2+10.1s+1)
%den = (s.^3+2*s.^2+10*s+9)

num = [1 10.1 1];
den = [1 2 10 9];

sys=tf(num, den);
%bode(sys)
%pzmap(sys)
%grid on


% Uppgift 2 b)

F=100;
Tmax=2.^13;
Ts=(2*pi)/F;
t=0:Ts:Tmax*Ts;

x1 = sin(t);
x2 = sin(3*t);
x3 = sin(5*t);

% Uppgift 2 c)

% Amplitud
y1 = lsim(sys,x1,t);
y2 = lsim(sys,x2,t);
y3 = lsim(sys,x3,t);

y1e = abs(evalfr(sys,1j))*x1;
y2e = abs(evalfr(sys,3j))*x2;
y3e = abs(evalfr(sys,5j))*x3;

max(y1)
% max(y1e)
% max(y2)
% max(y2e)
% max(y3)
% max(y3e)

% Fas
phi1 = angle(evalfr(sys,1j));
x1p = sin(t-phi1);
y1p = abs(evalfr(sys,1j))*x1p;

phi2 = angle(evalfr(sys,3j));
x2p = sin(t-phi2);
y2p = lsim(sys,x2p,t);

phi3 = angle(evalfr(sys,5j));
x3p = sin(t-phi3);
y3p = lsim(sys,x3p,t);

% Plot
subplot(2,1,1)
%plot(t, a1, 'r', t, a2, 'b');
plot(t,y_1)
hold on
plot(t,y1a)
axis([100 120 -5 5])
subplot(2,1,2);
plot(t, p, 'k');
axis([0 10 0 10]);



%% Uppgift 3 a)
clf
clc

t = -2*pi:0.01:2*pi;
y = square(t);
%plot(t,y, 'r', 'linewidth', 2)

% x(t)!=x(-t) => x är udda

% Uppgift 3 b)
for n = 1:5
    An = 0;
    Bn = 4*mod(n,2)/(n*pi);
    if(Bn ~= 0)
        disp(Bn)
    end
end

% Uppgift 3 c)
Fs=100;
k=0:1:2.^13-1;
wk=(2*pi*Fs*k)/8192;
ffy=fft(y, 8192);
plot(wk, abs(ffy))
hold on
grid on
