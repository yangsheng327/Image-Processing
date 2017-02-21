function [flamespeed, markstein]= KBLforderi(radius,speed,guess)



results = nlinfit(radius,speed, @KBLeqnderi, guess);

flamespeed = results(1);
markstein = results(2);

