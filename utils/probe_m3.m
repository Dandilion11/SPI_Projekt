function x = probe_m3(x,Vp)
i = 2:5;

x(i) = x(i).*(1 - Vp/x(1));
x(1) = x(1) - Vp;
end