M = 6;
Msp = 3;
Mr = 2*M;

tau = Msp/(M*M);

h = 2*pi / Mr;

x = 6.5*h; y = 2.5*h; z = 2*h;

m1 = 6; m2 = 2; m3 = 2;

x1 = m1*h; y1 = m2*h; z1 = m3*h;


E1 = exp( -( (x-x1)^2+(y-y1)^2+(z-z1)^2 )/(4*tau)  )

E2x = exp(pi*(x-x1)/(Mr*tau))

E2y = exp(pi*(y-y1)/(Mr*tau))

E2z = exp(pi*(z-z1)/(Mr*tau))

E2xl = E2x.^((-Msp+1) : Msp)
E2yl = E2y.^((-Msp+1) : Msp)
E2zl = E2z.^((-Msp+1) : Msp)

E3 = exp( -(pi*(0:Msp)/Mr).^2/tau)

E4 = exp(tau*(0:M/2).^2)

V0 = 1 * E1;
for k = -Msp+1 : Msp
   V1 = V0 * E2zl(k+Msp) * E3(abs(k)+1); 
   for j = -Msp+1 : Msp
      V2 = V1 * E2yl(j+Msp) * E3(abs(j)+1); 
      for i = -Msp+1 : Msp
          V3 = V2 * E2xl(i+Msp) * E3(abs(i)+1);
          
          
          a(i+Msp,j+Msp,k+Msp) = V3;
         
          b(i+Msp,j+Msp,k+Msp) = exp( -( (x-x1-h*i)^2 + (y-y1-h*j)^2+(z-z1-h*k)^2 )/(4*tau)  );
      end
   
   end
   
end

for k = -1
   V1 = V0 * E2zl(k+Msp) * E3(abs(k)+1); 
   for j = 0
      V1 = V1 * E2yl(j+Msp) * E3(abs(j)+1); 
      for i = -1
          V1 * E2xl(i+Msp) * E3(abs(i)+1)
      
       
         
        
      end
   
   end
   
end



