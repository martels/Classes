function rcfilter=rcfilter(t,data,a);
%
%
% inputs t = time axis as a vector
%        data = the data to be processed
%        a = the fraction of new signal to accept (dt/RC)
%
% output rcfilter = the filtered output, assuming it starts at
% zero. 
%
rcfilter=zeros(size(t));
rcfilter(1)=a*data(1);
for n=2:length(t);
  rcfilter(n)=(1-a)*rcfilter(n-1)+a*data(n);
end;
plot(t,data,t,rcfilter);grid on;
xlabel('t, time, sec');
ylabel('signals');
