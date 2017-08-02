function rhs=burgers_rhs6(t,u0t,dummy,k,eps)

u0=ifft(u0t); u0x=ifft(i*k.*u0t); u0xx=ifft(-(k.^2).*u0t); 

% model 6
% u_t = -(1.06351184)u*u_{x} + (0.57099496) u*u_{xx} - (0.59272620) u^2 * u_{xx}
rhs=  fft(-(1.06351184)*(u0x.*u0) +(0.57099496)*(u0xx.*u0)-(0.59272620)*(u0xx.*(u0.^2)) );


