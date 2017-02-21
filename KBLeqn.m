function time=KBLeqn(beta,radius)

flamespeed=beta(1);
markstein=beta(2);
constant=beta(3);

time=zeros(length(radius),1);
%t=zeros(length(r),1);
for i=1:length(radius);
    time(i,1)=(radius(i)+2*markstein*log(radius(i))-4*markstein^2/radius(i)-(8/3)*markstein^3*(radius(i)^(-2))-constant)/flamespeed;
end

