function N = A2N(A)

A = upper(A);

N = [];

letters = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

for i=1:length(A)
    
    J = 0;
for j=1:length(letters)
     
    if strcmp(char(letters(j)),A(i))
        J = j;
    end
    
end

  N(1,i) = J;

end


return
