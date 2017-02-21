function A = N2A(N)

if N<1
    N=1;
elseif N>26
    N=26;
end

letters = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};


A = char(letters(N));