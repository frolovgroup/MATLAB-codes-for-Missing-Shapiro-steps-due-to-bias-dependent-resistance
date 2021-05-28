clear all
b=0:0.01:1.16; %i_rf
z=length(b);
c=0:0.01:1.3; %i_dc
x=length(c);

cnt1=0;
cnt2=0;

h=6.626*10.^-34; %Planck's constant
f=2.63*10.^9; %rf frequency in Hz
e=1.6*10.^-19; %electron charge
Ic=3.3*10.^-6; %critical current
R=33; %Constant resistance
tau=(2*e*R*Ic)/(h);

Resistance=zeros(z,x);

% Code for choosing i_dc for the odd steps
for n=1:105

    avgphidot=zeros(1,x);
    phi_initial=0;
    phidot_initial=0;
    
    for m=1:x
        
        step=0.0003*10.^-9;
        t=0:0.0003*10.^-9:4.8*10.^-9; %t in second
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
    
    % Choosing i_dc at which change is slope occurs in I-V curves
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
    
    % Calculating I_bias for odd steps
    p=1;
    for i=1:5
        c_odd(i)=(c_steps(p)+c_steps(p+1))/2;
        p=p+2;
    end
    
    Res=@(Ibias)(33-(18*(0.03.^2)/((0.03.^2)+((Ibias-c_odd(1)).^2)))-(7.5*(0.03.^2)/((0.03.^2)+((Ibias-c_odd(2)).^2)))-(6.5*(0.03.^2)/((0.03.^2)+((Ibias-c_odd(3)).^2)))-(4*(0.03.^2)/((0.03.^2)+((Ibias-c_odd(4)).^2)))-(3*(0.03.^2)/((0.03.^2)+((Ibias-c_odd(5)).^2))));
 
    for m=1:x
        Resistance(n,m)= Res(c(m));
    end
    
    cnt1=cnt1+1
end

%Manually defining R(i_dc) for i_rf>1.04
%For i_rf=1.05 to 1.07
Res1=@(Ibias)(33-(18*(0.03.^2)/((0.03.^2)+((Ibias-0.05).^2)))-(7.5*(0.03.^2)/((0.03.^2)+((Ibias-0.24).^2)))-(6.5*(0.03.^2)/((0.03.^2)+((Ibias-0.41).^2)))-(4*(0.03.^2)/((0.03.^2)+((Ibias-0.575).^2)))-(3*(0.03.^2)/((0.03.^2)+((Ibias-0.72).^2))));
 
for n=106:108 
    for m=1:x
        Resistance(n,m)= Res1(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.08 to 1.09
Res2=@(Ibias)(33-(18*(0.03.^2)/((0.03.^2)+((Ibias-0.05).^2)))-(7.5*(0.03.^2)/((0.03.^2)+((Ibias-0.21).^2)))-(6.5*(0.03.^2)/((0.03.^2)+((Ibias-0.385).^2)))-(4*(0.03.^2)/((0.03.^2)+((Ibias-0.55).^2)))-(3*(0.03.^2)/((0.03.^2)+((Ibias-0.715).^2))));
               
for n=109:110 
    for m=1:x
        Resistance(n,m)= Res2(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.1 to 1.12
Res3=@(Ibias)(33-(18*(0.03.^2)/((0.03.^2)+((Ibias-0.05).^2)))-(7.5*(0.03.^2)/((0.03.^2)+((Ibias-0.18).^2)))-(6.5*(0.03.^2)/((0.03.^2)+((Ibias-0.36).^2)))-(4*(0.03.^2)/((0.03.^2)+((Ibias-0.525).^2)))-(3*(0.03.^2)/((0.03.^2)+((Ibias-0.685).^2))));
 
for n=111:113 
    for m=1:x
        Resistance(n,m)= Res3(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.13 to 1.14
Res4=@(Ibias)(33-(18*(0.03.^2)/((0.03.^2)+((Ibias-0.05).^2)))-(7.5*(0.03.^2)/((0.03.^2)+((Ibias-0.16).^2)))-(6.5*(0.03.^2)/((0.03.^2)+((Ibias-0.34).^2)))-(4*(0.03.^2)/((0.03.^2)+((Ibias-0.505).^2)))-(3*(0.03.^2)/((0.03.^2)+((Ibias-0.665).^2))));
 
for n=114:115 
    for m=1:x
        Resistance(n,m)= Res4(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.15 to 1.16
Res5=@(Ibias)(33-(18*(0.03.^2)/((0.03.^2)+((Ibias-0.05).^2)))-(7.5*(0.03.^2)/((0.03.^2)+((Ibias-0.15).^2)))-(6.5*(0.03.^2)/((0.03.^2)+((Ibias-0.315).^2)))-(4*(0.03.^2)/((0.03.^2)+((Ibias-0.485).^2)))-(3*(0.03.^2)/((0.03.^2)+((Ibias-0.645).^2))));
               
for n=116:117 
    for m=1:x
        Resistance(n,m)= Res5(c(m));
    end
    cnt1=cnt1+1
end

% Applying R(i_dc)
r=1;
Vy=zeros(1,z*x);
xax=zeros(1,z*x);

for n=1:z
    
    phi_initial_new=0;
    phidot_initial_new=0;
    
    for m=1:x
        
        xax(1,r)=b(n);
        
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
        
        Vy(1,r)=mean(phidotnew)/(f);
        phi_initial_new=phi(k);
        phidot_initial_new=phidotnew(k);
        r=r+1;
    end
    
    cnt2=cnt2+1
end

his = [transpose(xax),transpose(Vy)];
hist3(his,'Ctrs',{0:0.01:1.2 -0.2:0.2:20.2},'CDataMode','auto','EdgeColor','interp');
colormap('hot')
xlabel('$i_{rf}$', 'interpreter', 'latex');
ylabel('$V/[hf/2e]$', 'interpreter', 'latex');
yticks([0 2 4 6 8 10 12 14 16 18 20]);
cl=colorbar
lim = caxis
caxis([0 30])
cl.Title.String = "Bin Counts";
view(2)
xlim([0,1.16])
ylim([0,10.2])