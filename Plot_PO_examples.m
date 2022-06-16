orbit = load('ORBIT.DAT');

%figure
plot(orbit(:,1), orbit(:,2), 'LineWidth', 2, 'Color', 'black');
hold on
scatter(orbit(1,1), orbit(1,2),'filled','red','LineWidth',3);
%[x, y, z] = ellipsoid(0,0,0,0.867986935857794,0.846287262461350,0.824587589064904,30);
%surf(x, y, z)
set(gca,'xlim',[min(orbit(:,1)) max(orbit(:,1))]);
set(gca,'ylim',[min(orbit(:,2)) max(orbit(:,2))]);
xlabel('X');
ylabel('Y');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)