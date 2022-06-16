%close all;
orbit = load('ORBIT.DAT');

% figure
% plot3(orbit(:,1), orbit(:,2),orbit(:,3), 'LineWidth', 2, 'Color', 'black');
% hold on
% scatter3(orbit(1,1), orbit(1,2),orbit(1,3),'filled','red','LineWidth',3);
% %[x, y, z] = ellipsoid(0,0,0,1,0.975,0.95,30);
% %surf(x, y, z)
% xlabel('X');
% ylabel('Y');
% box('on');
% grid on
% set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

figure
plot(orbit(:,1), orbit(:,2), 'LineWidth', 2, 'Color', 'black');
hold on
%scatter(orbit(1,1), orbit(1,2),'filled','red','LineWidth',3);
[x, y, z] = ellipsoid(0,0,0,0.867986935857794,0.846287262461350,0.824587589064904,100);
surf(x, y, z)
set(gca,'xlim',[min(orbit(:,1)) max(orbit(:,1))]);
set(gca,'ylim',[min(orbit(:,2)) max(orbit(:,2))]);
xlabel('X');
ylabel('Y');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)


% figure;
% subplot(3,1,1);
% plot(orbit(:,7),orbit(:,1),'color','black','LineWidth',2)
% title('X');
% 
% subplot(3,1,2);
% plot(orbit(:,7),orbit(:,2),'color','black','LineWidth',2)
% title('Y');
% 
% subplot(3,1,3);
% plot(orbit(:,7),orbit(:,3),'color','black','LineWidth',2)
% title('Z');

% orbitN = load('Force.DAT');
% figure
% plot3(orbitN(:,1), orbitN(:,2),orbitN(:,3), 'LineWidth', 2, 'Color', 'black');
% %hold on
% %scatter3(orbit(1,1), orbit(1,2),orbit(1,3),'filled','red','LineWidth',3);
% %[x, y, z] = ellipsoid(0,0,0,1,0.975,0.95,30);
% %surf(x, y, z)
% xlabel('X');
% ylabel('Y');
% box('on');
% grid on
% set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

