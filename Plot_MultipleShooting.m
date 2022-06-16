orbitN = load('ORBIT_Num.DAT');
orbitA = load('ORBIT_Ana.DAT');
orbitR = load('ORBIT_Ref.DAT');
orbitAN = load('ORBIT_AnaRef.DAT');

figure
plot3(orbitA(:,1), orbitA(:,2), orbitA(:,3), 'LineWidth', 2, 'Color', 'blue');
hold on
%plot3(orbitN(:,1), orbitN(:,2), orbitN(:,3), 'LineWidth', 2, 'Color', 'black');
%plot3(orbitAN(:,1), orbitAN(:,2), orbitAN(:,3), 'LineWidth', 2, 'Color', 'green');
scatter3(orbitR(:,1), orbitR(:,2), orbitR(:,3),'filled','red','LineWidth',3);
% X label
xlabel('x');
% Y label
ylabel('y');
% Z label
zlabel('z');
grid on
set(gca,'FontName','Helvetica','FontSize',30)
%[x, y, z] = ellipsoid(0,0,0,0.766880337068647,0.751843467714360,0.677336457400324,50);
%surfl(x, y, z);

figure
%plot(orbitN(:,1), orbitN(:,2), 'LineWidth', 2, 'Color', 'black');

plot(orbitA(:,1), orbitA(:,2), 'LineWidth', 2, 'Color', 'blue');
hold on
%plot(orbitAN(:,1), orbitAN(:,2), 'LineWidth', 2, 'Color', 'green');
scatter(orbitR(:,1), orbitR(:,2),'filled','red','LineWidth',3);
% X label
xlabel('x');
% Y label
ylabel('y');
grid on
set(gca,'FontName','Helvetica','FontSize',30)