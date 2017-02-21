function [fittime, fitradius,instantspeed]=globalpolynomial(time,radius,polyorder,offset)

fittime=[];
fitradius=[];
instantspeed=[];

p=polyfit(time((1+offset):(length(time)-offset)),radius((1+offset):(length(radius)-offset)),polyorder); %Determines coefficients

pprime=polyderivative(p);

for i=1:length(time)
    xvec=[];
    x=time(i);
    for j=1:polyorder+1;
        xvec(j)=x^(polyorder-j+1); %[x^n x^(n-1)...x 1]
    end
    instantspeed(i)=pprime*xvec';
    fitradius(i)=p*xvec';
    fittime(i)=time(i);
end