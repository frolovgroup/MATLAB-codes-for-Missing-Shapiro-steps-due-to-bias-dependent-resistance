clear all
b=0.9;%i_rf
c=0:0.01:1.3; %i_dc
x=length(c);

h=6.626*10.^-34;%Planck's constant
f=2.63*10.^9;%rf frequency in Hz
e=1.6*10.^-19;%electron charge
Ic=3.3*10.^-6;%critical current
R=33;%constant (background) resistance
tau=(2*e*R*Ic)/(h);

% Code for choosing i_dc for the odd steps
avgphidot=zeros(1,x);%average voltage
phi_initial=0;
phidot_initial=0;
    
for m=1:x
        
    step=0.0003*10.^-9;
    t=0:0.0003*10.^-9:4.8*10.^-9; %t in seconds
    k=length(t);
        
    phi=zeros(1,k);
    phi(1)=phi_initial;
    phidotnew=zeros(1,k);
    phidotnew(1)=phidot_initial;
      
    %RK4 method
    for i=1:k-1
            
        phidot=@(t,phi)(tau*(b*sin(f*t)-sin(phi)+c(m)));
            
        k1 = phidot(t(i),phi(i));
        k2 = phidot(t(i)+0.5*step,phi(i)+0.5*k1*step);
        k3 = phidot(t(i)+0.5*step,phi(i)+0.5*k2*step);
        k4 = phidot(t(i)+step,phi(i)+k3*step);

        phi(i+1) = phi(i)+((k1+2*k2+2*k3+k4)/6)*step;
        phidotnew(i+1)=phidot(t(i+1),phi(i+1));
            
    end
        
    avgphidot(1,m)=mean(phidotnew)/(f);
    phi_initial=phi(k);
    phidot_initial=phidotnew(k);
end
   
% Choosing i_dc at which slope changes in I-V curve
dVdI_choose=zeros(1,x);
c_steps=zeros(1,10);
p=1;

for i=1:x-1
    dVdI_choose(i)=(avgphidot(i+1)-avgphidot(i))/(c(i+1)-c(i));
end
    
for i=1:x-1
    if i==1
        if ceil(dVdI_choose(i+1)-dVdI_choose(i))>40
            c_steps(p)=c(i+1);
            p=p+1;
        end
    elseif i==x-1
        if ceil(dVdI_choose(i)-dVdI_choose(i-1))>40
            c_steps(p)=c(i+1);
            p=p+1;
        end
    else
        if ceil(dVdI_choose(i+1)-dVdI_choose(i))>40 && ceil(dVdI_choose(i)-dVdI_choose(i-1))<40
            c_steps(p)=c(i+1);
            p=p+1;
        end
    end
end
    
% Calculating i_dc for odd steps
c_odd=zeros(1,5);
p=1;
for i=1:5
    c_odd(i)=(c_steps(p)+c_steps(p+1))/2;
    p=p+2;
end
  
%Applying R(i_dc)
V=zeros(1,x);%average voltage
phi_initial_new=0;
phidot_initial_new=0;
 
for m=1:x
        
    step=0.0003*10.^-9;
    t=0:0.0003*10.^-9:4.8*10.^-9; %t in seconds
    k=length(t);
        
    phi=zeros(1,k);
    phi(1)=phi_initial_new;
    phidotnew=zeros(1,k);
    phidotnew(1)=phidot_initial_new;
        
    %Resistance peak at Shapiro step 1 with HWHM=0.04 and Peak height=22.5 ohms
    Res=@(Ibias)(33+((22.5*(0.04.^2)/((0.04.^2)+((Ibias-c_odd(1)).^2)))));
        
    for i=1:k-1
            
        phidot=@(t,phi)(tau*(b*sin(f*t)-sin(phi)+c(m)));
            
        k1 = phidot(t(i),phi(i));
        k2 = phidot(t(i)+0.5*step,phi(i)+0.5*k1*step);
        k3 = phidot(t(i)+0.5*step,phi(i)+0.5*k2*step);
        k4 = phidot(t(i)+step,phi(i)+k3*step);

        phi(i+1) = phi(i)+((k1+2*k2+2*k3+k4)/6)*step;
        phidotnew(i+1)=phidot(t(i+1),phi(i+1));
        tau=(2*e*Res(c(m))*Ic)/(h);
            
    end
        
    V(1,m)=mean(phidotnew)/(f);
    phi_initial_new=phi(k);
    phidot_initial_new=phidotnew(k);
end

plot(V(1,1:x),c)
xlabel('$V/[hf/2e]$', 'interpreter', 'latex');
ylabel('$i_{dc}$', 'interpreter', 'latex')