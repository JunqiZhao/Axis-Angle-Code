function[x] = denominator(Bx, By, Bz, Gx, Gy, Gz)
x1 = Bx * cos(asin(-Gx/((Gx.^2+Gy.^2+Gz.^2).^(1/2))));
x2 = By * sin(atan(-Gy/Gz)) * sin(asin(-Gx/((Gx.^2+Gy.^2+Gz.^2).^(1/2))));
x3 = Bz * cos(atan(-Gy/Gz)) * sin(asin(-Gx/((Gx.^2+Gy.^2+Gz.^2).^(1/2))));
x = x1+x2+x3;
end

