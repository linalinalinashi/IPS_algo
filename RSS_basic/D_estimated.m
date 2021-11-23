function [De_a,De_b,De_c]  = D_estimated(Adet,m,pra,prb,prc,P_LED,h)

global P_LED

%calculate 3 distances with 3 LED
De_a = power(((m+1)*Adet*(h)^(m+1)*P_LED(1))/(2*pi*pra), 1/(m+3));
De_b = power(((m+1)*Adet*(h)^(m+1)*P_LED(2))/(2*pi*prb), 1/(m+3));
De_c = power(((m+1)*Adet*(h)^(m+1)*P_LED(3))/(2*pi*prc), 1/(m+3));

end

