%% script for EECE 2150 Lab 13 Fall 2016
%%
%% Original created by N, McGruer, March, 2016
%% Modified by D. Brooks, Nov 2016, Ali Vakili and N. McGruer Nov. 2017
%%
%% This script creates a few periods of a square wave, finds a Fourier
%% representation of that signal using the Matlab "fft" function, and then
%% uses samples in frequency of a 1st order low-pass analog (CT) filter to
%% simulate filtering the signal by appropriately modifying the FFT
%% coefficients, then using an inverse FFT to compute the filtered signal
clc
clear all
close all
% First make a time vector to simulate time samples to get 10 period of 
% the signal, so the time vector is made based on frequency of the signal, 
% with samping frequency of 10kHz

f_sig_in = 100;

Fs = 100*f_sig_in;
t=[0:1/Fs:10/f_sig_in];

L = length(t);

% Second make a square wave with the same length, with the given frequency
% so in case of f = 10, there will be 10 periods  in one second.

% Here we simply create a sinusoid with amplitude 0.9 added to a DC constant
% of 1 and then "quantize" it to either 1 or 0 using the Matlab "floor"
% function.
% 
% There are many other ways to create a square wave, for example by using
% the Matlab "zeros" and "ones" functions and manipulating vector indices

% Note that "floor" makes any value between 1 and 2 equal to 1 and any
% %value between 0 and 1 equal to 0, resulting is a square
% wave between 0 and 1.

Sig_in= floor(1 + 0.9*sin(2*pi*f_sig_in*t));

% Plot the square wave.

figure(1)
plot(t,Sig_in,'linewidth',3)
axis([0 t(end) -0.2, 1.2])
title('Square Wave, Time Domain')
xlabel('Time, seconds')
ylabel('Amplitude (volts, for example)')
grid on
pause

% Compute a Fourier representation of the square wave using the Matlab 
% implementation of the Discrete Fourier Transform in its "fft: function. 
% We need to scale the results after fft function to get the correct
% amplitude

fft_sig_in=fft(Sig_in);

fft_sig_in=fft_sig_in/L;

% Plot the magnitude of the elements of the fft, as computed by MATLAB.  
% The FFT returns complex values in general, so we use the "abs" function
% to plot the magnitudes.
%
% For reasons that are too complicated to describe here but certainly worth
% a discussion in class or during the lab, low frequencies are at both ends
% of the frequency range and high frequencies in the middle!

% To plot, we first construct a vector of frequencies (f) to match the
% elements of the fft.  

f = [0:1:L-1]*Fs/L;
 
stem(f,abs(fft_sig_in),'linewidth',3), grid on
xlabel('Frequency, Hz')
title('Amplitude of Fourier Representation of Square Wave Segment')
pause

%Now just look at the lowest frequencies to see more detail

axis([0 0.1*Fs 0 1.1*max(abs(fft_sig_in))])
pause

% Now to be sure we know where we are actually computing values, we will mark
% each point with an 'x;

%hold on
%plot(f,abs(fft_sig_in),'.r','markersize',20), grid on
%pause

% Electrical engineers are generally used to looking at frequency with DC
% in the center and higher positive / negative frequencies going out to the
% right and left. However by convention the DFT is computed with DC as the
% first coefficient, and the periodicity of DT frequencies causes those low
% frequencies so show up at the end of the DFT as well. 

% To get around this, Matlab has a function "fftshift" to center the FFT
% around DC. Here we will use it to plot the magnitude of the elements of
% the fft, as we normally look at them.   

% Again, we construct a vector of frequencies to match the elements of the
% fft after the shift. Note that as described, zero is at the center, and will be
% associated with the zero %frequency, or DC component of the Fourier
% representation 
 
hold off
f2 = [-(L-1)/2:1:(L)/2]*Fs/L;
stem(f2,abs(fftshift(fft_sig_in)),'linewidth',3), grid on
xlabel('Frequency, Hz')
title('Amplitude of Fourier Representation of Square Wave Segment')
pause

% Now again just look at the lowest frequencies.

axis([-.1*Fs 0.1*Fs 0 1.1*max(abs(fftshift(fft_sig_in)))])
pause

% Now we use the inverse transform implemented in the "ifft" function to go
% back to the time domain and plot to demonstrate that the fft function
% followed by the ifft function gives the original time domain waveform
% back.

fft_sig_scaled = fft_sig_in*L;

sig_in_ifft=ifft(fft_sig_scaled);
plot(t,sig_in_ifft,'r','linewidth',3), grid on
xlabel('Time, seconds')
title('Amplitude of Reconstructed Signal in Red')
axis([0 t(end) -.2 1.2])
hold on
plot(t,Sig_in,'--b', 'linewidth',3), grid on
title('Amplitude of Reconstructed Signal in Red and Original in Blue ')
legend('Reconstructed Signal','Original Signal')
pause

% Now we will filter the frequency domain version with frequency samples of
% a low-pass filter:

% H(w)=-1/(1+jwRC)

% Here, let's assume that R=1e5 and C=1e-7, so H(w)=-1/(1+j0.01w)
% It's a little tricky doing the filtering this way because the fft has
% two sides, so this will look a little complicated, but the idea is to
% treat both sides the same, knowing that one is the complex conjugate
% of the other, and that the first term is the DC or a0 term.  
% As we go on we will do more with MATLAB and learn how all this works...

% To keep DC in the center we use the f2 vector to supply the frequencies
% to use. So the following statement creates a vector Hw containing the
% values of H(w) computed at the frequencies in f2.

hold off
Hw=-2./(1+1i*.01*2*pi.*f2);

% Plot the transfer function, also called frequency response, of the
% filter, using the axis command to plot only the part near the center of
% the frequency spectrum. 

% Note that H(w) is also complex-valued so we need to use the abs() function
% again to plot its magnitude

plot(f2, abs(Hw),'.-','linewidth',3), grid on
axis([-.1*Fs .1*Fs 0 max(abs(Hw))])
xlabel('Frequency, Hz')
title('Amplitude of Samples of Low Pass Transfer Function')

pause

% Now filter the square wave, using the conventional (not Matlab) view of 
% the fft of the square wave.  Note that the filter transfer function is
% set up to operate on frequency components with zero in the middle, so
% we have to do it this way to match the fft to Hw.

fft_sig_in_filtered=fftshift(fft_sig_in).*Hw;

% Now plot the frequency domain representation of the 
% filtered and unfiltered square waves.
% first, the unfiltered version, axes adjusted to see the center of the 
% frequency spectrum (the part near zero frequency).

stem(f2,abs(fftshift(fft_sig_in)),'linewidth',3), grid on
axis([-.1*Fs .1*Fs 0 max(abs(fftshift(fft_sig_in)))])
xlabel('Frequency, Hz')
title('Amplitude of Fourier Representation of Original Square Wave Segment')
hold on
%Now the filtered version

stem(f2,abs(fft_sig_in_filtered),'--r','linewidth',3), grid on
xlabel('Frequency, Hz')
title('Amplitude of Filtered  Fourier Representation of Square Wave Segment')
legend('Original Signal','Filtered Signal')
pause
hold off

% Finally, convert the frequency domain representation back to the
% time domain and plot. Note we have to "un-fftshift" to get DC back to the
% beginning and end so that it matches the way the FFT function is set up.

% first the filtered version

fft_sig_in_filtered_scaled = fft_sig_in_filtered*L;

%Sig_in_filter=-1*ifft(ifftshift(fft_sig_in_filtered_scaled));
Sig_in_filter=ifft(ifftshift(fft_sig_in_filtered_scaled));
plot(t,Sig_in_filter,'linewidth',3), grid on
axis([0 t(end) -1.2 1.2])
xlabel('Time, seconds')
title('Amplitude of Filtered Signal')
hold on
% Now, the original time-domain waveform for comparison

plot(t,Sig_in,'r','linewidth',3), grid on
xlabel('Time, seconds')
title('Amplitude of Original and Filtered Signals')
legend('Filtered Signal','Original Signal')
hold off
