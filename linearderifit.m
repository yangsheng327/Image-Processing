function [flamespeed,markstein]= linearderifit(sb, kappa)

p = polyfit(kappa,sb,1);

    flamespeed=p(2);
    markstein=-p(1);

