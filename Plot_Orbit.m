orbitN = load('Orbit1.DAT');
EP = orbitN(1,:);
orbitN(1,:) =[];

figure
scatter3(orbitN(1,1), orbitN(1,2), orbitN(1,3),'filled','green','LineWidth',3);
hold on
plot3(orbitN(:,1), orbitN(:,2), orbitN(:,3), 'LineWidth', 2, 'Color', 'black');
scatter3(EP(1),EP(2),EP(3),'filled','red','LineWidth',3)
xlabel('x');
ylabel('y');
zlabel('z');
grid on
set(gca,'FontSize',15,'FontWeight','bold','LineWidth',1)
[x, y, z] = ellipsoid(0,0,0,0.766880337068647,0.751843467714360,0.677336457400324,50);
surfl(x, y, z);

figure
scatter(orbitN(1,1), orbitN(1,2),'filled','green','LineWidth',3);
hold on
plot(orbitN(:,1), orbitN(:,2), 'LineWidth', 2, 'Color', 'black');
scatter(EP(1),EP(2),'filled','red','p')
xlabel('x','FontWeight','bold');
ylabel('y','FontWeight','bold');
grid on
box on
set(gca,'FontSize',30,'FontWeight','bold','LineWidth',2)