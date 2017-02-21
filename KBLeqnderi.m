function speed = KBLeqnderi(beta,radius)

flamespeed=beta(1);
markstein=beta(2);

%speed=zeros(length(radius),1);
speed=zeros(length(radius),1);
for i=1:length(radius);
    speed(i) = flamespeed./( 1+ 2*markstein./radius(i) + 4*markstein^2./radius(i).^2 + (16/3)*markstein^3./radius(i).^3);
end


