%% Uppgift 3.1 a)

% An = 0
% Bn = 4/(k*pi)
% Se fk.tex

% Uppgift 3.1 b)
clf
clc
syms n x k

T=1;
w=2*pi/T;
M=200;
x=0;
t=T*(0:M-1)/M;

for n=1:100
    Ak = 0;
    Bk = 4*mod(n,2)/(n*pi);
    x = x + Ak*cos(n*w*t) + Bk*sin(n*w*t);
end
plot(t,x)

%% Uppgift 3.2 a)
clf
clc

% y(t) = g(w)sin(wt+phi(w))
% num  = (s^2+10.1s+1)
% den  = (s^3+2*s^2+10*s+9)

num = [1 10.1 1];
den = [1 2 10 9];
Gs=tf(num, den);

% Uppgift 3.2 b)
F=100;
N=2.^13;
Ts=1/F;
Tmax = (N-1)*Ts;
t=0:Ts:Tmax; 

% Insignals
x1 = sin(t);
x2 = sin(3*t);
x3 = sin(5*t);

% % Plots
% subplot(3,1,1)
% plot(t, x1)
% legend('sin(t)')
% axis([0 80 -1 1])
% subplot(3,1,2)
% plot(t, x2)
% legend('sin(3t)')
% axis([0 80 -1 1])
% subplot(3,1,3)
% plot(t, x3)
% legend('sin(5t)')
% axis([0 80 -1 1])

% Uppgift 3.2 c)

% y = x sent through sys using lsim
y1 = lsim(Gs,x1,t);
y2 = lsim(Gs,x2,t);
y3 = lsim(Gs,x3,t);

% y = correct amplitude, wrong phase
y1e = abs(evalfr(Gs,1j))*x1;
y2e = abs(evalfr(Gs,3j))*x2;
y3e = abs(evalfr(Gs,5j))*x3;

% y = correct amp and phase due to eq 2
phi1 = angle(evalfr(Gs,1j));
x1p = sin(t+phi1);
y1p = abs(evalfr(Gs,1j))*x1p;

phi2 = angle(evalfr(Gs,3j));
x2p = sin(3*t+phi2);
y2p = abs(evalfr(Gs,3j))*x2p;

phi3 = angle(evalfr(Gs,5j));
x3p = sin(5*t+phi3);
y3p = abs(evalfr(Gs,5j))*x3p;

% % % Plots
% % a)
% bode(sys)
% pzmap(sys)

% c)
subplot(3,1,1)
plot(t,y1,'--b', t,y1p,':r', t, x1,'k')
axis([0 30 -1 1.2])
legend('lsim1', 'ekv 2', 'x1'), title('w=1')

subplot(3,1,2)
plot(t,y2, '--b', t, y2p, ':r',t, x2,'k')
axis([0 30 -5 5])
legend('lsim2', 'ekv 2','x2'), title('w=3')

subplot(3,1,3)
plot(t,y3, '--b', t, y3p, ':r',t,x3,'k')
axis([0 30 -1.3 1.2])
legend('lsim3', 'ekv 2','x3'), title('w=5')

%% Uppgift 3.3 a)
clf
clc

N = 2.^13;
F = 100;
Ts = 1/F;
t = 0:Ts:40*pi;
x = square(t);
% x(t) != x(-t) => x är udda

% Uppgift 3.3 b)
fprintf('3.3 b\n')
fprintf('FK enligt ekv 1:\n\n')
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
kf=@(wk) (N*wk)/(2*pi*F);
ffx=fft(x, N);

% Uppgift 3.3 d)
B1 = (2*abs(ffx(k+1)))/N;
fprintf('3.3 d\n')
fprintf('FK enligt fft (ekv 10):\n\n')
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

%plot(t, y);

% Enligt ekv 8
fprintf('3.3 e\n')
fprintf('FK enligt ekv 8:\n\n')
disp(abs(evalfr(H,1j))*4/(pi*1))
disp(abs(evalfr(H,3j))*4/(pi*3))
disp(abs(evalfr(H,5j))*4/(pi*5))

% Enligt fft
ffy = fft(y, N);
B2 = (2*abs(ffy(k+1)))/N;

fprintf('FK enligt fft (ekv 10):\n\n')
B2max1 = max((B2(1:ceil(kf(2)))));
B2max2 = max((B2(ceil(kf(2)):ceil(kf(4)))));
B2max3 = max((B2(ceil(kf(4)):ceil(kf(6)))));
disp(B2max1)
disp(B2max2)
disp(B2max3)

% Plots
subplot(3,1,1)
plot(t, x, 'r', 'linewidth', 1.5)
axis([-0 10 -1.5 1.5])
grid on
subplot(3,1,2)
plot(wk, abs(ffx))
axis([0 6 0 6000])
grid on
subplot(3,1,3)
plot(wk, abs(B2))
axis([0 6 -0 1.5])
grid on

%% Uppgift 3.4 a)
clc
clf

syms s;

Ts = s*(s-1j)*(s+1j)*(s-5j)*(s+5j)*(s-7j)*(s+7j)*(s-9j)*(s+9j);
num = [1 0 156 0 7374 0 106444 0 99225 0];
Tp = tf(num, 1);
% % Plots
%bode(Tp);
%grid on

% Uppgift 3.4 b: n=10, c: n=11)
Np = 1;
for n = 1:12
    Np = Np*(s+4);
end
den = sym2poly(Np);
sys = tf(num,den);
w={1,5000};
% % Plots
%bode(sys,w);
%grid on

% Uppgift 3.4 d)
scale = abs(evalfr(sys, 3j));
sys2 = tf(num/scale, den);
abs(evalfr(sys2,3j));
% % Plots
% bode(sys2, 'r');
% grid on
% hold on
% bode(sys, '--b');

% Uppgift 3.4 e)
N=8192;
k = 0:(N-1);
wk = (2*pi*F*k)/(N);
F = 100;
Ts = 1/F;
t = 0:Ts:(N-1)*Ts;
x = square(t);

% x=square(t)
yx = lsim(sys2,x,t);
ffy = fft(yx, 8192);
By = (2*abs(ffy(k+1)))/N;
% % Plots
% plot(t,x, 'k', t, yx, 'b');
% legend('x(t)', 'y(t)')
% axis([0 30 -1.5 1.5])
% plot(wk, abs(By));
% axis([0 20 0 .5])


% % Amplitud hos Notchfiltrets utsignal
% fprintf('Notch-amp square:')
% disp(max(yx))
% % Amplitud genom DTF och FFT
% fprintf('FFT-amp square:')
% disp(max(By))

% ysignalen
num = [1 10.1 1];
den = [1 2 10 9];
H=tf(num, den);
y = lsim(H,x,t);
yy = lsim(sys2, y,t);
ffy = fft(yy, 8192);
By = (2*abs(ffy(k+1)))/N;
% % Plots
% plot(t,y, 'k', t, yy);
% legend('x(t)', 'y(t)')
% axis([0 30 -3 3])
plot(wk, abs(By));
axis([0 20 0 1.5])

% % Amplitud hos Notchfiltrets utsignal
% fprintf('Notch-amp yy:')
% disp(max(yy))
% % Amplitud genom DTF och FFT
% fprintf('FFT-amp yy:')
% disp(max(By))