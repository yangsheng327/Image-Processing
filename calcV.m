function V=calcV(radius, beta)
%Function calcV calculates the derivative matrix
% v_{np}=\frac{\partial f(x_n,theta)} {\partial \theta_p} \right|_{theta^o}

rmin=min(radius);
for i=1:length(radius)
    
    tau=findtau(radius(i)/(exp(beta(2))-rmin/exp(1)));
    
    V(i,1)=(exp(beta(2))-rmin/exp(1))*exp(beta(3));
    V(i,2)=(exp(beta(3))*exp(beta(2)))*((1/tau^2*log(tau))-real(expint(log(tau^2)))+beta(1)+((rmin/exp(1)-exp(beta(2)))/(radius(i)*tau^3*(log(tau))^2)));
    V(i,3)=(exp(beta(2))-rmin/exp(1))*((1/(tau^2*log(tau)))-real(expint(log(tau^2)))+beta(1))*exp(beta(3));

end
