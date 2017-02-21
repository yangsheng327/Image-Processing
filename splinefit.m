function [fittime, fitradius,instantspeed] = splinefit(time,radius)

fittime=time;
fitradius=radius;
instantspeed=[];


% delta=0.000001;

% pp = csaps(radius,time,0.999); %smoothing spline

%axes(handles.axes1)
%plot(radius,ppval(pp,radius),'k-')


% for i=1:length(radius)
%     instantspeed(i)=delta/(ppval(pp,radius(i)+delta)-ppval(pp,radius(i)));
% end


delta=0.000000000001;

pp = csaps(time,radius,0.99999999999); %smoothing spline

for i=1:length(time)
     instantspeed(i)=(ppval(pp,time(i)+delta)-ppval(pp,time(i)))/delta;
 end