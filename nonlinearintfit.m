function [flamespeed, markstein, constant] = nonlinearintfit(radius,time, initial)

   if isnan(initial(1))
   initial(1) = 200;
   end
   if isnan(initial(2))
   initial(2) = 0.01;
   end
   if isnan(initial(3))
   initial(3) = 0.5;
   end

   
   
rmin=min(radius);
tolerance=1E-7;

%Initial Guess for beta matrix
beta(1,1)=(initial(1)*initial(3))/(-2*initial(2));
if((rmin/exp(1)-2*initial(2))>0)
    beta(2,1)=log(rmin/exp(1)-2*initial(2));
else
    beta(2,1)=-2;
end
beta(3,1)=log(1/initial(1));  


    constant=exp(beta(1));
    flamespeed=1/exp(beta(3));
    markstein=(exp(beta(2))-radius(1)/exp(1))/(-2);

    
delta(1,1)=10;
delta(2,1)=10;
delta(3,1)=10; %Initialize step so we can begin iterations



format long e;
k=0;
lambda=2; %Initialize lambda (start at 2 since we are going to halve the guess)


% h=figure;


while( max(abs((lambda*delta)./beta)) > tolerance )
    V=calcV(radius,beta); %Calculates the derivative matrix

    [Q R]=qr(V); %Splits the derivative matrix into an orthogonal and a diagonal matrix
    R1=R(1:length(beta),1:length(beta));
    Q1=Q(:,1:length(beta));
    
    for i=1:length(radius)
        tau=findtau(radius(i)/(exp(beta(2))-rmin/exp(1)));
        tpred(i)=exp(beta(3))*((1/(tau^2*log(tau)))-real(expint(log(tau^2)))+beta(1))*(exp(beta(2))-rmin/exp(1));
        z(i,1)=time(i)-tpred(i); %Residual vector
    end
    
    S0=sum(z.*z);

    w=Q1'*z;

    rsize = size(R1);

    for ri=1:rsize(1)
        for rj=1:rsize(2)
            if isnan(R1(ri,rj))
                R1(ri,rj)=0;
            end
        end
    end
    
    delta=pinv(R1)*w;
    
%     D=eye(3);
%     D(1,1)=V(1,1);
%     D(2,2)=V(2,2);
%     D(3,3)=V(3,3);
%     delta=(V'*V+k*eye(3))\(V'*z);

    lambda=2;
    
    S1=S0+1;
    beta1=[0; 0; 0];
    while(S1>S0)
        lambda=lambda/2;
        beta1=beta+lambda*delta;
        for i=1:length(radius)
            tau=findtau(radius(i)/(exp(beta1(2))-rmin/exp(1)));
            tpred(i)=exp(beta1(3))*((1/(tau^2*log(tau)))-real(expint(log(tau^2)))+beta1(1))*(exp(beta1(2))-rmin/exp(1));
            z1(i,1)=time(i)-tpred(i); %Residual vector
        end 
        S1=sum(z1.*z1);
    end
    
    
% 
%     hold on
%     cla
%     plot(radius,time,'b.')
%     plot(radius,tpred,'r-')
%     drawnow;
%     
    
    beta=beta1;
    
    constant=exp(beta1(1));
    flamespeed=1/exp(beta1(3));
    markstein=(exp(beta1(2))-radius(1)/exp(1))/(-2);
    
   %[constant flamespeed markstein]

end

%Convert from tilde form to nontilde

beta(2)=exp(beta(2))-radius(1)/exp(1);
beta(3)=exp(beta(3));

constant=beta(1);
flamespeed=1/beta(3);
markstein=beta(2)/(-2);







