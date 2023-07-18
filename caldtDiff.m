function dt_diffusion = caldtDiff(T,numSegs,b_norm,a)

    k = 8.617e-5;  % eV * K-1
    kt = k*T;   %eV
    R = (sqrt(3)+2)/4*numSegs*b_norm;
    n_sia = pi*sqrt(3)*R^2/a^2;
    Em = 0.06+0.07*n_sia^(1/3);  %eV
    D0_original = 0.009*n_sia^(-0.621); %cm^2 * s^-1
    D0 = D0_original*1e16;  % A^2 * s^-1
    Dn = D0*exp(-Em/kt);     % A^2 * s^-1
    dt_diffusion = b_norm^2/(2*Dn);   % s

end

