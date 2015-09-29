%% Uppgift 1 a)

T=1;
w=2*pi/T;
M=200;
t=T*(0:M-1)/M;
y=sin(w*t);
plot(t,y)

% Svar:
%      An=1/T*(int(cos(w*n*x), x, 0, T/2) - int(cos(w*n*x), x, T/2, T));
%      Bn=1/T*(int(sin(w*n*x), x, 0, T/2) - int(sin(w*n*x), x, T/2, T));

%% Uppgift 1 b)
clf
syms n x k

T=1;
w=2*pi/T;
M=200;
y=0;
f=0;
t=T*(0:M-1)/M;

for n=1:10
    An=1/T*(int(cos(w*n*x), x, 0, T/2) - int(cos(w*n*x), x, T/2, T));
    Bn=1/T*(int(sin(w*n*x), x, 0, T/2) - int(sin(w*n*x), x, T/2, T));
    y=y + An*cos(n*w*t) + Bn*sin(n*w*t);
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
bode(sys)
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

%g = evalfr(sys, 2)
y1 = lsim(sys,x1,t)
y2 = lsim(sys,x2,t);
y3 = lsim(sys,x3,t);

plot(y1);
asdf = freqresp(sys, 1)*x1;
max(asdf)
%plot(t,y)

%% Uppgift 3 a)
clf
clc

t = -2*pi:0.01:2*pi;
y = square(t);
plot(t,y,'linewidth', 2)

% x(t)!=x(-t) => x är udda

% Uppgift 3 b)
% FS(x)= 
w=1;
T=2*pi;
n=0:5;
    An=1/T*(int(-cos(w*n*x), x, -pi, 0) + int(cos(w*n*x), x, 0, pi));
    Bn=1/T*(int(-sin(w*n*x), x, -pi, 0) + int(sin(w*n*x), x, 0, pi));

% Uppgift 3 c)
Fs=1000;
k=0:1:2.^13;
wk=(2*pi)./(Fs*k);
ffy=fft(y,2^nextpow2(2.^13))/2.^13;
plot(ffy)
hold on
grid on
