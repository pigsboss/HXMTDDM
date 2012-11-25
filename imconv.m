function g=imconv(f,h)
g=real(ifft2(fft2(ifftshift(h)).*fft2(f)));
return