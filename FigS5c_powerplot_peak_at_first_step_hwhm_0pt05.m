clear all
b=0:0.01:1.2; %i_rf
z=length(b);
c=0:0.01:1.3; %I_dc
x=length(c);

dVdI=zeros(z,x); %dV/dI

cnt1=0;
cnt2=0;

h=6.626*10.^-34; %Planck's constant
f=2.63*10.^9; %rf frequency in Hz
e=1.6*10.^-19; %electron charge
Ic=3.3*10.^-6; %critical current
R=33; %Constant resistance
tau=(2*e*R*Ic)/(h);

Resistance=zeros(z,x);%R(i_dc) for different i_rf
       
for n=1:105

    % Code for choosing I_bias for the odd steps
    avgphidot=zeros(1,x); 
    phi_initial=0;
    phidot_initial=0;
    
    %RK4 method for solving non-linear differential equation
    for m=1:x
        
        step=0.0003*10.^-9;
        t=0:0.0003*10.^-9:4.8*10.^-9; %t in seconds
        k=length(t);
        
        phi=zeros(1,k);
        phi(1)=phi_initial;
        phidotnew=zeros(1,k);
        phidotnew(1)=phidot_initial;
        
        for i=1:k-1
            
            phidot=@(t,phi)(tau*(b(n)*sin(f*t)-sin(phi)+c(m)));
            
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
    c_odd=zeros(1,5);
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
    p=1;
    for i=1:5
        c_odd(i)=(c_steps(p)+c_steps(p+1))/2;
        p=p+2;
    end
    
    Res=@(Ibias)(33+((22.5*(0.05.^2)/((0.05.^2)+((Ibias-c_odd(1)).^2)))));
    for m=1:x
        Resistance(n,m)= Res(c(m));
    end
    
    cnt1=cnt1+1
end

%Manually defining R(i_dc) for i_rf>1.04
%For i_rf=1.05 to 1.2
Res1=@(Ibias)(33+((22.5*(0.05.^2)/((0.05.^2)+((Ibias-0.05).^2)))));

for n=106:121 
    for m=1:x
        Resistance(n,m)= Res1(c(m));
    end
    cnt1=cnt1+1
end

% Applying R(i_dc)
    
for n=1:z
    
    V=zeros(1,x); %Average voltage 
    phi_initial_new=0;
    phidot_initial_new=0;
    
    for m=1:x
        
        step=0.0003*10.^-9;
        t=0:0.0003*10.^-9:4.8*10.^-9; %t in second
        k=length(t);
        
        phi=zeros(1,k);
        phi(1)=phi_initial_new;
        phidotnew=zeros(1,k);
        phidotnew(1)=phidot_initial_new;
        
        for i=1:k-1
            
            phidot=@(t,phi)(tau*(b(n)*sin(f*t)-sin(phi)+c(m)));
            
            k1 = phidot(t(i),phi(i));
            k2 = phidot(t(i)+0.5*step,phi(i)+0.5*k1*step);
            k3 = phidot(t(i)+0.5*step,phi(i)+0.5*k2*step);
            k4 = phidot(t(i)+step,phi(i)+k3*step);

            phi(i+1) = phi(i)+((k1+2*k2+2*k3+k4)/6)*step;
            phidotnew(i+1)=phidot(t(i+1),phi(i+1));
            
            tau=(2*e*Resistance(n,m)*Ic)/(h);
            
        end
        
        V(1,m)=mean(phidotnew)/(f);% Average V/[hf/2e]
        phi_initial_new=phi(k);
        phidot_initial_new=phidotnew(k);
    end
    
    %Calculation of differential resistance
    for q=1:x
        if q==1
            dVdI(n,q)=((V(q+1)-V(q))/(c(q+1)-c(q)));
        elseif q==x
            dVdI(n,q)=((V(q)-V(q-1))/(c(q)-c(q-1)));
        else
            dVdI(n,q)=((V(q+1)-V(q-1))/(c(q+1)-c(q-1)));
        end
    end
    
    cnt2=cnt2+1
    
end

[B,C]=meshgrid(b,c);
pcolor(B,C,transpose(dVdI));
colormap('hot')
shading interp 
cl=colorbar
lim = caxis
caxis([0 100])
cl.Title.String = "dV/dI";
xlabel('$i_{rf}$', 'interpreter', 'latex');
ylabel('$i_{dc}$', 'interpreter', 'latex');
xlim([0,1.2])