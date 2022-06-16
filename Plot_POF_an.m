%load('POF.mat')
%POT = POT1_an;

POT1 = POT_an1;%POT(find(POT(:,8)<127.737157294961),:);
POT2 = POT_an2;%POT(find(POT(:,8)>=127.737157294961),:);
POT = [POT_an2;POT_an1];
%POT = [flipud(POT_an2);POT_an1];

figure
subplot(2,1,1)
plot(POT1(:,8),POT1(:,2), 'LineWidth', 2, 'Color', 'black');
hold on
plot(POT1(:,8),POT1(:,14), 'LineWidth', 2, 'Color', 'blue');
plot(POT1(:,8),POT1(:,15), 'LineWidth', 2, 'Color', 'blue');
pic01 = fill([POT1(:,8);flipud(POT1(:,8))],[POT1(:,15);flipud(POT1(:,14))],'b');
set(pic01,'edgealpha', 0, 'facealpha', 0.4);

plot(POT2(:,8),POT2(:,2), 'LineWidth', 2, 'Color', 'black');
plot(POT2(:,8),POT2(:,14), 'LineWidth', 2, 'Color', 'blue');
plot(POT2(:,8),POT2(:,15), 'LineWidth', 2, 'Color', 'blue');
pic02 = fill([POT2(:,8);flipud(POT2(:,8))],[POT2(:,15);flipud(POT2(:,14))],'b');
set(pic02,'edgealpha', 0, 'facealpha', 0.4);
ylabel('y_p');
box('on');
grid on;
set(gca,'xticklabel',[]);
set(gca,'xlim',[20 200]);
set(gca,'ylim',[0.8 1.25]);
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2,'position',[0.15 0.55 0.75 0.4])

subplot(2,1,2)
plot(POT(:,8),POT(:,12), 'LineWidth', 2, 'Color', 'blue','DisplayName','\omega_s');
hold on
plot(POT(:,8),POT(:,11), 'LineWidth', 2, 'Color', 'black','DisplayName','\omega_l');
plot(POT(:,8),POT(:,13), 'LineWidth', 2, 'Color', 'red','DisplayName','\omega_v');
plot(POT(:,8),2.*pi./POT(:,7), 'LineWidth', 2, 'Color', 'green','DisplayName','\omega');
xlabel('Tn');
ylabel('\omega');
box('on');
grid on;
set(gca,'xlim',[20 200]);
set(gca,'ylim',[0.7 1]);
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2,'position',[0.15 0.15 0.75 0.4])