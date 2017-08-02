function rhs=burgers_rhs5(t,u0t,dummy,k,eps)

u0=ifft(u0t); u0x=ifft(i*k.*u0t); u0xx=ifft(-(k.^2).*u0t); u0xxx=ifft(-i*(k.^3).*u0t);

% model 5
% u_t = -(1.02242009) u* u_{x} + (0.08725293) u_{xx} + (0.09602962)* u * u_{xx} - (0.09726498) u^2 *u_{xx}
rhs= (0.08725293)*(-k.^2).*u0t ...
+ fft(-(1.02242009)*(u0x.*u0) + (0.09602962)*(u0xx.*u0) -(0.09726498)*(u0xx.*(u0.^2)) );


