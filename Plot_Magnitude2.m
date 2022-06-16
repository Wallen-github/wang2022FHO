
load('G:\Master\科研工作\SpaceCraftInAsteroid\数据统计\AsteroidWithSatellite\user_data.mat');
[VMBA,validM] = data_pro(MBA);
[VMCA,valid] = data_pro(MCA);
[VNEA,valid] = data_pro(NEA);
[VTNA,valid] = data_pro(TNA);
[VJTA,valid] = data_pro(JTA);
total = [MBA;MCA;NEA;TNA;JTA];
%total(find(isnan(total(:,6))),:)=[];
%total(find(isnan(total(:,7))),:)=[];
[Vtotal,total] = data_pro(total);
total = [total,Vtotal];


% figure;
% scatter(VMBA(:,1),VMBA(:,2),'filled','b','o')
% hold on
% scatter(VMCA(:,1),VMCA(:,2),'filled','r','s')
% scatter(VNEA(:,1),VNEA(:,2),'filled','m','d')
% scatter(VTNA(:,1),VTNA(:,2),'filled','k','^')
% scatter(VJTA(:,1),VJTA(:,2),'filled','c','p')
% xlabel('$\overline{V}^{3rd}$','interpreter','latex');
% ylabel('$\overline{V}^{SRP}$','interpreter','latex');
% %box('on');
% grid on
% set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2)
% 
% axes('position',[0.42 0.42 0.45 0.45]);
% scatter(VMBA(:,1),VMBA(:,2),'filled','b','o')
% hold on
% scatter(VMCA(:,1),VMCA(:,2),'filled','r','s')
% scatter(VNEA(:,1),VNEA(:,2),'filled','m','d')
% scatter(VTNA(:,1),VTNA(:,2),'filled','k','^')
% scatter(VJTA(:,1),VJTA(:,2),'filled','c','p')
% x= linspace(0,1.5E-3,100);
% y = x;
% y_low = 0.01*x;
% y_high = 100*x;
% plot(x,y,'color','black','LineWidth',1)
% plot(x,y_low,'color','blue','LineWidth',1,'LineStyle','--')
% plot(x,y_high,'color','blue','LineWidth',1,'LineStyle','--')
% xlim([0 1.5E-3])
% ylim([0 1.5E-3])
% box('on');
% grid on
% set(gca,'FontSize',15,'FontWeight','bold','LineWidth',1)

figure;
x= linspace(min(min(Vtotal)),max(max(Vtotal)),100);
y = x;
y_low = 0.1*x;
y_high = 10*x;
scatter(VMBA(:,1),VMBA(:,2),'filled','b','o','DisplayName','Main belt members')
hold on
scatter(VMCA(:,1),VMCA(:,2),'filled','r','s','DisplayName','Mars crossers')
scatter(VNEA(:,1),VNEA(:,2),'filled','m','d','DisplayName','Near Earth objects')
%scatter(VTNA(:,1),VTNA(:,2),'filled','k','^','DisplayName','Trans-Neptunian objects')
scatter(VJTA(:,1),VJTA(:,2),'filled','c','p','DisplayName','Jupiter Trojans')
plot(x,y,'color','black','LineWidth',1,'DisplayName','\gamma = 1')
plot(x,y_low,'color','blue','LineWidth',1,'LineStyle','--','DisplayName','\gamma_{min} = 0.1')
plot(x,y_high,'color','blue','LineWidth',1,'LineStyle','--','DisplayName','\gamma_{max} = 10')
xlabel('$\overline{V}^{3rd}$','interpreter','latex');
ylabel('$\overline{V}^{SRP}$','interpreter','latex');
xlim([min(Vtotal(:,1)) max(Vtotal(:,1))])
ylim([min(Vtotal(:,2)) max(Vtotal(:,2))])
box('on');
%grid on
set(gca,'FontSize',18,'FontWeight','bold','LineWidth',2,'YScale','log','XScale','log')

load('example_in_paper.mat');
[Vforpaper,data] = data_pro(dataforpaper);

scatter(Vforpaper(:,1),Vforpaper(:,2),'filled','b','p','DisplayName','examples')

function [V,data] = data_pro(data)
data(find(isnan(data(:,1))),:)=[];
data(find(isnan(data(:,5))),:)=[];
data(find(isnan(data(:,4))),:)=[];
data(find(isnan(data(:,10))),:)=[];
data(find(isnan(data(:,11))),:)=[];
data(find(isnan(data(:,12))),:)=[];

GG = 6.67E-11;
AU = 1.496E11;
SM = 0.01;
kappa = 1.44;
PSRP = 4.5605E-6;
mS = 2E30;

rA = data(:,4)/2*1000;
rB = data(:,5)/2*1000;
rho = data(:,1)*1000;
TA = data(:,12)*60*60;
mA = 4.0/3.*pi.*rA.^3.*rho;
mB = 4.0/3.*pi.*rB.^3.*rho;
r=(GG.*mA.*(TA./(2.*pi)).^2).^(1.0./3);
rS = data(:,10)*AU;
r1 = data(:,11)*1000;
rhoS = kappa.*SM.*PSRP.*rS.^2;
V3rd = GG.*mB.*r.^2./r1.^3;
VSRP = rhoS.*r./rS.^2;
V0 = GG.*mA./r;

V = [V3rd./V0,VSRP./V0];
end


