[Y, FS]=audioread('Recording.wav');

y = Y(1:90000);
p = 1/numel(y);

t = 0:p:1000;
t = t(1:90000);


v1 = fft(y);

figure; plot(t, abs(v1), '-k');


figure;plot(t*1e3,v1,'b-');grid on;
title('Pulse');
xlabel('t, Time, msec');ylabel('v, Volts');

f=fftaxisshift(fftaxis(t));
V1=fftshift(fft(v1));

figure;cxplot(f/1e3,V1);grid on;
title('Pulse');
xlabel('f, Freq, kHz');ylabel('v, Volts/Hz');
%{
% zoom to center
figure;cxplot(f/1e3,V1);grid on;
xzoom(-10,10);
title('Pulse');
xlabel('f, Freq, kHz');ylabel('v, Volts/Hz');
%}
