restart
f0 = d^2-c1^2-c2^2-c3^2;
f1 = l1*l2*l3*c3^2*a^2*b^2+k^3;
f2 = c3^2*a^2*b^2*(l1*l2+l1*l3+l2*l3)-k^2*(d^2-a^2-b^2);
f3 = c3^2*a^2*b^2*(l1+l2+l3) - k*(d^2*a^2+c1^2*b^2-c1^2*a^2+c3^2*b^2-a^2*b^2);
I = ideal(f0,f1,f2,f3)
end
