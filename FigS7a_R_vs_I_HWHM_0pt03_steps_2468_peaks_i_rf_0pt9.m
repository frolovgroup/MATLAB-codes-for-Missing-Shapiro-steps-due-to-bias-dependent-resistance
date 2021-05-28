%Code to plot R(I) at i_rf=0.9 for resistance peaks at steps 2,4,6,8
clear all
ibias=0:0.001:1.1; %i_dc
s=length(ibias);
R=zeros(1,s);%R(i_dc)

Res=@(ibias)(33+(11.5*(0.03.^2)/((0.03.^2)+((ibias-0.28).^2)))+(7*(0.03.^2)/((0.03.^2)+((ibias-0.445).^2)))+(5*(0.03.^2)/((0.03.^2)+((ibias-0.605).^2)))+(4*(0.03.^2)/((0.03.^2)+((ibias-0.75).^2))));

for i=1:s
    R(i)=Res(ibias(i));
end
plot(ibias,R)
xlabel('$i_{dc}$', 'interpreter', 'latex');
ylabel('$R (\Omega)$', 'interpreter', 'latex');
xlim([0,1])
ylim([0,60])