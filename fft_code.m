clc;
clear all;
close all;
f1 = input('Enter first frequency');
f2 = input('Enter second frequency');
T = 1/3000;
t = 0:T:1-T;
x = sin(2*pi*f1*t);
M=input('Enter the value of N');
%N indicates the number of times the dft is split to compute fft (it
%decides accuracy of fft generated)
z = fft_split(x,M);

subplot(2,1,1);
plot(t,x);
xlabel('time');
ylabel('amplitude');
title('Input Signal');

subplot(2,1,2);
plot(abs(z));
axis([0 820 0 1000]);
xlabel('frequency');
ylabel('magnitude');
title('Magnitude Plot');

