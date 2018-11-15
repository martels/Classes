function h=cxplot(t,v);
% Plot a complex function - 
% real in blue -o,
% imag in  red-x,
% abs in green -
% by Chuck DiMarzio
% Northeastern University
% November 2018
%
% input t for horizontal axis, v for vertical
% output h is a handle to the plot.
% To bring the plot to the foreground enter figure(h);
% Then use functions like xlabel, ylabel, etc. 

h=plot(t,real(v),'b-o',t,imag(v),'r-x',t,abs(v),'g-');
legend('Real','Imag','Abs');