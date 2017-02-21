function tau=findtau(r)

%Iterate to find tau given r
%r=1/(tau*ln(tau))

e=exp(1);
tol=1E-5; %Error tolerance
h=0.001; %Step size in calculating the derivative
error=1; %Initialize error to high value so that it is greater than tolerance and begins loops

if(r<0 && r>-e)
    1
    display('error - this value of r is not allowed')
    tau=1/e;
    tauguess=0;
elseif r<0
    tauguess=0.5;
    while error>tol
        y1=r*tauguess*log(tauguess)-1;
        y2=r*(tauguess+h)*log(tauguess+h)-1;
        taunew=tauguess-(y1*h)/(y2-y1);
        if(taunew<(1/e))
            taunew=1/e;
        elseif(taunew>1)
            taunew=1-h-tol;
        end
        error=abs(taunew-tauguess);
        tauguess=taunew;
    end
    
else
    tauguess=2;
    while error>tol
        y1=r*tauguess*log(tauguess)-1;
        y2=r*(tauguess+h)*log(tauguess+h)-1;
        taunew=tauguess-(y1*h)/(y2-y1);
        if(taunew<1)
            taunew=1;
        end
        error=abs(taunew-tauguess);
        tauguess=taunew;
    end
end
tau=tauguess;
