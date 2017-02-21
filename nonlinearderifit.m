function [flamespeed, markstein]= nonlinearderifit(s, kappa, guess)

if isnan(guess(1))
    guess(1) = 200;
end
if isnan(guess(2))
    guess(2) = 0; 
end


markstein=0;
flamespeed=0;


results = nlinfit(s,kappa, @nlinfunction, guess);


flamespeed=results(1);
markstein=results(2);


