function limits=xzoom(x1,x2);
% xzoom(x1,x2) to zoom plot in x direction only between x2 and x2
% by Chuck DiMarzio, November 2018, Northeastern University
%

moose=axis;
axis([x1,x2,moose(3:4)]);
limits=axis;
