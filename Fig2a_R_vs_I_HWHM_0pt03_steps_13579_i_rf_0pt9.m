%Code to plot R(i_dc) at i_rf=0.9
clear all
ibias=0:0.001:1.1; %i_dc
s=length(ibias);
R=zeros(1,s);%R(i_dc)

Res=@(ibias)(33+(21.5*(0.03.^2)./((0.03.^2)+((ibias-0.195).^2)))+(7.5*(0.03.^2)./((0.03.^2)+((ibias-0.37).^2)))+(6.5*(0.03.^2)./((0.03.^2)+((ibias-0.53).^2)))+(5*(0.03.^2)./((0.03.^2)+((ibias-0.685).^2)))+(4*(0.03.^2)./((0.03.^2)+((ibias-0.825).^2))));

for i=1:s
    R(i)=Res(ibias(i));
end
plot(ibias,R)
xlabel('$i_{dc}$', 'interpreter', 'latex');
ylabel('$R (\Omega)$', 'interpreter', 'latex');
xlim([0,1])
ylim([0,60])