function [flamespeed, markstein, constant]= KBLforint(radius,time,guess)


results = nlinfit(radius,time, @KBLeqn, guess);

flamespeed = results(1);
markstein = results(2);
constant = results(3);

