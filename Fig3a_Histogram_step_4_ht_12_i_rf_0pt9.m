clear all
b=0.9;%i_rf=0.9
ht=12;%Height of resistance peak at Shapiro step 4
c=0:0.01:1.3; %i_dc
x=length(c);

h=6.626*10.^-34;
f=2.63*10.^9;
e=1.6*10.^-19;
Ic=3.3*10.^-6;
R=33;
tau=(2*e*R*Ic)/(h);

hwhm=0:0.01:0.5;%HWHM of resistance peak at Shapiro step 4
s=length(hwhm);

% Code for choosing i_dc for the odd steps
avgphidot=zeros(1,x);
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
    
%Choosing i_dc at which change in slope occurs in I-V curve
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

% Applying R(i_dc)
Vy=zeros(1,s*x);
xax=zeros(1,s*x);
r=1;
cnt=0;

for j=1:s
    
    phi_initial_new=0;
    phidot_initial_new=0;
    
    for m=1:x
        
        xax(1,r)=hwhm(j);
        
        step=0.0003*10.^-9;
        t=0:0.0003*10.^-9:4.8*10.^-9; %t in second
        k=length(t);
        
        phi=zeros(1,k);
        phi(1)=phi_initial_new;
        phidotnew=zeros(1,k);
        phidotnew(1)=phidot_initial_new;
        
        Res=@(Ibias)(33+((ht)*((hwhm(j).^2)/((hwhm(j).^2)+((Ibias-c_even(2)).^2)))));
        
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
        
        Vy(1,r)=mean(phidotnew)/(f);
        phi_initial_new=phi(k);
        phidot_initial_new=phidotnew(k);
        r=r+1;
    end
    cnt=cnt+1
end 

his = [transpose(xax),transpose(Vy)];
hist3(his,'Ctrs',{0:0.01:0.5 -0.2:0.2:20.2},'CDataMode','auto','EdgeColor','interp');
xlabel('$i_{rf}$', 'interpreter', 'latex');
ylabel('$V/[hf/2e]$', 'interpreter', 'latex');
yticks([0 2 4 6 8 10 12 14 16 18 20]);
colormap('hot')
cl=colorbar
lim = caxis
caxis([0 30])
cl.Title.String = "Bin Counts";
view(2)
xlim([0,0.5])
ylim([0.2,14.2])