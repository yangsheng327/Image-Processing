function P = polyperimeter(x,y)


%   P = polyperimeter(x,y) returns the area and the centroid
%   coordinates of the polygon specified by the vertices in 
%   the vectors x and y.

P = 0;

if length(x) ~= length(y)

            return;
            display('x and y do not have same length');

else

            for i=1:length(x)
                if i==1
                     P = P + ((x(length(x))-x(1))^2+(y(length(y))-y(1))^2)^0.5;
                else
                      P = P + ((x(i)-x(i-1))^2+(y(i)-y(i-1))^2)^0.5;
                end
            end
 
end