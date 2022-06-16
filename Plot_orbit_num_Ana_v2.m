orbitN = load('ORBIT_Numerical.DAT');
orbitA = load('ORBIT_Analytical.DAT');
orbitAN = load('ORBIT_AnaNum.DAT');

EP = orbitN(1,:);
orbitN(1,:) =[];
orbitA(1,:) =[];
orbitAN(1,:) =[];

figure
scatter3(EP(1),EP(2),EP(3),'filled','red','LineWidth',3)
hold on
plot3(orbitN(:,1), orbitN(:,2), orbitN(:,3), 'LineWidth', 2, 'Color', 'black');
plot3(orbitAN(:,1), orbitAN(:,2), orbitAN(:,3), 'LineWidth', 2, 'Color', 'blue');
% X label
xlabel('X');
% Y label
ylabel('Y');
% Z label
zlabel('Z');
grid on
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',1)
% 
% dR = sqrt(sum((orbitN(:,1:6)-orbitA(:,1:6)).^2,2));
% figure
% plot(orbitN(:,7),dR,'color','black','LineWidth',3)
% % 创建 ylabel
% ylabel('$D_{AN}$','FontWeight','bold','interpreter','latex');
% % 创建 xlabel
% xlabel('$time$','FontWeight','bold','interpreter','latex');
% grid on
% set(gca,'FontSize',20,'FontWeight','bold','LineWidth',1)
% 
% dR = sqrt(sum((orbitAN(:,1:6)-orbitA(:,1:6)).^2,2));
% figure
% plot(orbitAN(:,7),dR,'color','black','LineWidth',3)
% % 创建 ylabel
% ylabel('$D_{ANA}$','FontWeight','bold','interpreter','latex');
% % 创建 xlabel
% xlabel('$time$','FontWeight','bold','interpreter','latex');
% grid on
% set(gca,'FontSize',20,'FontWeight','bold','LineWidth',1)

dR = sqrt(sum((orbitAN(:,1:6)-orbitN(:,1:6)).^2,2));
meandR = mean(dR);
figure
p3 = plot(orbitAN(:,7),dR,'color','black','LineWidth',3);
hold on
plot([min(orbitAN(:,7)) max(orbitAN(:,7))],[meandR meandR],'color','black','LineWidth',3,'LineStyle','--');
% 创建 ylabel
ylabel('state vector variation','FontWeight','bold','interpreter','tex');
% 创建 xlabel
xlabel('time','FontWeight','bold','interpreter','tex');
grid on
legend([p3],'3rd X P5')
set(gca,'FontSize',20,'FontWeight','bold','LineWidth',1)
legend([p1,p2,p3])