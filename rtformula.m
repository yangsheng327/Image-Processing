function r = rtformula(beta, t)

r0 = beta(1);
A = beta(2);
alpha = beta(3);

r = r0 + A*(t.^alpha);

