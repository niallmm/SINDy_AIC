function rhs=burgers_rhs7(t,u0t,dummy,k,eps)

u0=ifft(u0t); u0x=ifft(i*k.*u0t); u0xx=ifft(-(k.^2).*u0t); u0xxx=ifft(-i*(k.^3).*u0t);

% model 7
% u_t = -(1.00965530) u * u_{x} + (0.10296556) u_{xx}
rhs= (0.10296556)*(-k.^2).*u0t + fft(-(1.00965530)*(u0x.*u0));



