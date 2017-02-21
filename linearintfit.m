function [flamespeed,markstein,constant]= linearintfit(radius, time)

low=1;
high=length(radius);


    for i=low:high
        A(i-low+1,:)=[radius(i) -1 2*log(radius(i))];
        B(i-low+1,:)=[time(i)];
    end
    x=inv(A'*A)*A'*B;

    flamespeed=1/x(1);
    markstein=x(3)*flamespeed;
    constant=x(2)*flamespeed;
