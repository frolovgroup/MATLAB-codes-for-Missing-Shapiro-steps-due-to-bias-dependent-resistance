clear all
b=0:0.01:1.2; %i_rf
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
    
    %Choosing i_dc at which change in slope occurs in I-V curves
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
    
    % Calculating i_dc for even steps
    c_even=zeros(1,4);
    p1=2;
    for i=1:4
        c_even(i)=(c_steps(p1)+c_steps(p1+1))/2;
        p1=p1+2;
    end
    
    Res=@(Ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((Ibias-c_even(1)).^2)))+(7*(0.03.^2)/((0.03.^2)+((Ibias-c_even(2)).^2)))+(5*(0.03.^2)/((0.03.^2)+((Ibias-c_even(3)).^2)))+(4*(0.03.^2)/((0.03.^2)+((Ibias-c_even(4)).^2))));
        
    for m=1:x
        Resistance(n,m)= Res(c(m));
    end
    
    cnt1=cnt1+1
end

%Manually defining R(i_dc) for i_rf>1.04
%For i_rf=1.05 to 1.06
Res1=@(Ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((Ibias-0.145).^2)))+(7*(0.03.^2)/((0.03.^2)+((Ibias-0.325).^2)))+(5*(0.03.^2)/((0.03.^2)+((Ibias-0.495).^2)))+(4*(0.03.^2)/((0.03.^2)+((Ibias-0.65).^2))));
        
for n=106:107 
    for m=1:x
        Resistance(n,m)= Res1(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.07 to 1.08
Res2=@(Ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((Ibias-0.13).^2)))+(7*(0.03.^2)/((0.03.^2)+((Ibias-0.31).^2)))+(5*(0.03.^2)/((0.03.^2)+((Ibias-0.475).^2)))+(4*(0.03.^2)/((0.03.^2)+((Ibias-0.635).^2))));
        
for n=108:109 
    for m=1:x
        Resistance(n,m)= Res2(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.09 to 1.11
Res3=@(Ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((Ibias-0.11).^2)))+(7*(0.03.^2)/((0.03.^2)+((Ibias-0.29).^2)))+(5*(0.03.^2)/((0.03.^2)+((Ibias-0.465).^2)))+(4*(0.03.^2)/((0.03.^2)+((Ibias-0.62).^2))));
        
for n=110:112
    for m=1:x
        Resistance(n,m)= Res3(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.12 to 1.13
Res4=@(Ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((Ibias-0.1).^2)))+(7*(0.03.^2)/((0.03.^2)+((Ibias-0.27).^2)))+(5*(0.03.^2)/((0.03.^2)+((Ibias-0.44).^2)))+(4*(0.03.^2)/((0.03.^2)+((Ibias-0.6).^2))));
        
for n=113:114 
    for m=1:x
        Resistance(n,m)= Res4(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.14 to 1.15
Res5=@(Ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((Ibias-0.1).^2)))+(7*(0.03.^2)/((0.03.^2)+((Ibias-0.25).^2)))+(5*(0.03.^2)/((0.03.^2)+((Ibias-0.425).^2)))+(4*(0.03.^2)/((0.03.^2)+((Ibias-0.585).^2))));
        
for n=115:116 
    for m=1:x
        Resistance(n,m)= Res5(c(m));
    end
    cnt1=cnt1+1
end

%For i_rf=1.16 to 1.17
Res6=@(Ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((Ibias-0.1).^2)))+(7*(0.03.^2)/((0.03.^2)+((Ibias-0.23).^2)))+(5*(0.03.^2)/((0.03.^2)+((Ibias-0.41).^2)))+(4*(0.03.^2)/((0.03.^2)+((Ibias-0.575).^2))));
        
for n=117:118 
    for m=1:x
        Resistance(n,m)= Res6(c(m));
    end
    cnt1=cnt1+1
end
    
%For i_rf=1.18 to 1.2
Res7=@(Ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((Ibias-0.1).^2)))+(7*(0.03.^2)/((0.03.^2)+((Ibias-0.22).^2)))+(5*(0.03.^2)/((0.03.^2)+((Ibias-0.395).^2)))+(4*(0.03.^2)/((0.03.^2)+((Ibias-0.555).^2))));
               
for n=119:121
    for m=1:x
        Resistance(n,m)= Res7(c(m));
    end
    cnt1=cnt1+1
end
        
% Applying R(i_dc)

Vy=zeros(1,z*x);
xax=zeros(1,z*x);
r=1;

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
xlabel('$i_{rf}$', 'interpreter', 'latex');
ylabel('$V/[hf/2e]$', 'interpreter', 'latex');
yticks([0 2 4 6 8 10 12 14 16 18 20]);
cl=colorbar
lim = caxis
caxis([0 30])
cl.Title.String = "Bin Counts";
view(2)
xlim([0,1.2])
ylim([0.2,10.5])