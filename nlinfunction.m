function kappa=nlinfunction(b,s)
flamespeed=b(1);
markstein=b(2);

for i=1:length(s);
    kappa(i,1)=-(flamespeed./(2*markstein)).*((s(i)./flamespeed).^2).*log((s(i)./flamespeed).^2);
      
end

