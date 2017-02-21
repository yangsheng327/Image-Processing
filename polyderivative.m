function pprime=polyderivative(p)
n=length(p)-1;
derivmatrix=zeros(length(p));
for i=1:n
    derivmatrix(i+1,i)=length(p)-i;
end

pprime=(derivmatrix*p')';