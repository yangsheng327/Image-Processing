function [fittime,fitradius,fitspeed] = localpolynomial(time,radius,polyorder,points)

fittime=[];
fitradius=[];
fitspeed=[];

if(points-2*floor(points/2) == 0) %if points is even
    points=points+1; %Forces odd number of points to be used
end

extend=(points-1)/2; %Number of points in each direction around central point
 

for i=1:length(radius)-2*extend %Cycle though each point but don't look at end points since they aren't valid
    index=i+extend; %index of current center
    p=polyfit(time(index-extend:index+extend),radius(index-extend:index+extend),polyorder); %Determines coefficients y=p_1*x^n + p_2*x^(n-1)+ ... +  p_{n+1}
    

    pprime=polyderivative(p);
    xvec=[];
    x=time(index);
    for j=1:polyorder+1;
        xvec(j)=x^(polyorder-j+1); %[x^n x^(n-1)...x 1]
    end
    fitspeed(i)=pprime*xvec';
    fitradius(i)=p*xvec';
    fittime(i)=x;
    
   
end
