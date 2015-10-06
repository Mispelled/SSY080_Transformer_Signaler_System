%% Uppgift 3.1 a)

% An = 0
% Bn = 4/(k*pi)
% Se uppgift1.txt

% Uppgift 3.1 b)
clf
clc
syms n x k

T=1;
w=2*pi/T;
M=200;
x=0;
f=0;
t=T*(0:M-1)/M;

for n=1:100
    An = 0;
    Bn = 4*mod(n,2)/(n*pi); %mod(n,2) för att få jämna k.
    x = x + An*cos(n*w*t) + Bn*sin(n*w*t);
end
plot(t,x)
hold on

%% Uppgift 3.2 a)
clf
clc

%num = (s.^2+10.1s+1)
%den = (s.^3+2*s.^2+10*s+9)

num = [1 10.1 1];
den = [1 2 10 9];
sys=tf(num, den);

% Uppgift 3.2 b)

F=100;
N=2.^13;
Ts=1/F;
Tmax = (N-1)*Ts;
t=0:Ts:Tmax; 

x1 = sin(t);
x2 = sin(3*t);
x3 = sin(5*t);

% Uppgift 3.2 c)

% Amplitud
% Red
y1 = lsim(sys,x1,t);
y2 = lsim(sys,x2,t);
y3 = lsim(sys,x3,t);

% Black
y1e = abs(evalfr(sys,1j))*x1;
y2e = abs(evalfr(sys,3j))*x2;
y3e = abs(evalfr(sys,5j))*x3;

% Fas
% Green
phi1 = angle(evalfr(sys,1j));
x1p = sin(t+phi1);
y1p = abs(evalfr(sys,1j))*x1p;

phi2 = angle(evalfr(sys,3j));
x2p = sin(3*t+phi2);
y2p = abs(evalfr(sys,3j))*x2p;

phi3 = angle(evalfr(sys,5j));
x3p = sin(5*t+phi3);
y3p = abs(evalfr(sys,5j))*x3p;

% Plots

% a)
% bode(sys)
% pzmap(sys)

% c)
subplot(3,1,1)
plot(t,y1, ':r', t, y1e, 'k', t, y1p, '--b')
axis([0 30 -1 1.2])
legend('y1', 'y1e', 'y1p')

subplot(3,1,2)
plot(t,y2, ':r', t, y2e, 'k', t, y2p, '--b')
axis([0 30 -5 5])
legend('y2', 'y2e', 'y2p')

subplot(3,1,3)
plot(t,y3, ':r', t, y3e, 'k', t, y3p, '--b')
axis([0 30 -1.3 1.2])
legend('y3', 'y3e', 'y3p')

grid on
%% Uppgift 3.3 a)
clf
clc

N = 2.^13;
F = 100;
Ts = 1/F;
t = 0:Ts:40*pi; % Ändrade denna för att få en "hårdare plot". 
                % Förändrar det något för våra Ck?
x = square(t);
% x(t) != x(-t) => x är udda

% Uppgift 3.3 b)
fprintf('Fourierkoefficienter enligt ekv 1:\n\n')
for n = 1:5;
    An = 0;
    Bn = 4*mod(n,2)/(n*pi);
    if(Bn)
        disp(Bn)
    end
end

% Uppgift 3.3 c)
k = 0:(N-1);
wk = (2*pi*F*k)/(N);
ffx=fft(x, N);

% % Plots
%  subplot(2,1,1)
%  plot(t, x, 'r', 'linewidth', 1.5)
%  axis([-0 10 -1.5 1.5])
%  subplot(2,1,2)
%  plot(wk, abs(ffx))
%  axis([0 7 0 6000])
%  hold on
%  grid on

% Uppgift 3.3 d)
B1 = (2*abs(ffx(k+1)))/N;
fprintf('Fourierkoefficienter enligt fft:\n\n')
B1max1 = max((B1(1:ceil(kf(2)))));
B1max2 = max((B1(ceil(kf(2)):ceil(kf(4)))));
B1max3 = max((B1(ceil(kf(4)):ceil(kf(6)))));
disp(B1max1)
disp(B1max2)
disp(B1max3)


% Uppgift 3.3 e)
num = [1 10.1 1];
den = [1 2 10 9];
H=tf(num, den);

y = lsim(H,x,t);

% Enligt ekv 8
fprintf('FK enligt ekv 8:\n\n')
disp(abs(evalfr(H,1j))*4/(pi*1))
disp(abs(evalfr(H,3j))*4/(pi*3))
disp(abs(evalfr(H,5j))*4/(pi*5))

% Enligt fft
ffy = fft(y, N);
B2 = (2*abs(ffy(k+1)))/N;
kf=@(wk) (N*wk)/(2*pi*F);

fprintf('Enligt fft:\n\n')
B2max1 = max((B2(1:ceil(kf(2)))));
B2max2 = max((B2(ceil(kf(2)):ceil(kf(4)))));
B2max3 = max((B2(ceil(kf(4)):ceil(kf(6)))));
disp(B2max1)
disp(B2max2)
disp(B2max3)
 
% Plot
subplot(2,1,1)
plot(wk, abs(B1))
axis([0 6 -0 1.5])
grid on
subplot(2,1,2)
plot(wk, abs(B2))
axis([0 6 -0 1.5])
grid on