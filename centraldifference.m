function [fittime,fitradius,instantspeed]=centraldifference(time,radius)

fitradius=[];
instantspeed=[];

for i=1:length(radius)-2
    index=i+1;
    instantspeed(i)=(radius(index+1)-radius(index-1))/(time(index+1)-time(index-1));
    fitradius(i)=radius(index);
    fittime(i)=time(index);
end
