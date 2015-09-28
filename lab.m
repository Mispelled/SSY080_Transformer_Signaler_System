%% Uppgift 1 a)

T=1;
w=2*pi/T;
M=200;
t=T*(0:M-1)/M;
y=sin(w*t);
plot(t,y)

%Svar: an = 1/T_0(integral_0^T/2(sin(nx)dx) + integral_T/2^T(-sin(nx)dx)
%      bn = 1/T_0(integral_0^T/2(cos(nx)dx) + integral_T/2^T(-cos(nx)dx)

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
    plot(t,y)
    hold on
end