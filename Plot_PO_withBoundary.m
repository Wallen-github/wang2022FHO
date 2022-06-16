POT = load('POF.DAT');
%load('POF_withBon.mat');
%POT = POT_re_withBon;

figure
subplot(2,1,1)
hold on
%POT = POT_an_withBon1;
plot(POT(:,8),POT(:,14), 'LineWidth', 2, 'Color', 'blue');
plot(POT(:,8),POT(:,15), 'LineWidth', 2, 'Color', 'blue');
pic01 = fill([POT(:,8);flipud(POT(:,8))],[POT(:,15);flipud(POT(:,14))],'b');
set(pic01,'edgealpha', 0, 'facealpha', 0.4);
plot(POT(:,8),POT(:,2), 'LineWidth', 2, 'Color', 'black');

% hold on
% POT = POT_an_withBon2;
% plot(POT(:,8),POT(:,14), 'LineWidth', 2, 'Color', 'blue');
% plot(POT(:,8),POT(:,15), 'LineWidth', 2, 'Color', 'blue');
% pic01 = fill([POT(:,8);flipud(POT(:,8))],[POT(:,15);flipud(POT(:,14))],'b');
% set(pic01,'edgealpha', 0, 'facealpha', 0.4);
% plot(POT(:,8),POT(:,2), 'LineWidth', 2, 'Color', 'black');
% POT = [POT_an_withBon2;POT_an_withBon1];

xlabel('Tn');
ylabel('y_p');
box('on');
grid on
set(gca,'xticklabel',[]);
set(gca,'xlim',[min(POT(:,8)) max(POT(:,8))]);
set(gca,'ylim',[0.8 1.2]);
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2,'position',[0.15 0.55 0.75 0.4])

subplot(2,1,2)
plot(POT(:,8),POT(:,12), 'LineWidth', 2, 'Color', 'blue','DisplayName','\omega_s');
hold on
%plot(POT(:,8),POT(:,11), 'LineWidth', 2, 'Color', 'black','DisplayName','\omega_l');
%plot(POT(:,8),POT(:,13), 'LineWidth', 2, 'Color', 'red','DisplayName','\omega_v');
plot(POT(:,8),2.*pi./POT(:,7), 'LineWidth', 2, 'Color', 'green','DisplayName','\omega');
xlabel('Tn');
ylabel('\omega');
box('on');
grid on;
set(gca,'xlim',[min(POT(:,8))  max(POT(:,8))]);
set(gca,'ylim',[0.8 1.6]);
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2,'position',[0.15 0.15 0.75 0.4])

figure
plot(POT(:,9),POT(:,12), 'LineWidth', 2, 'Color', 'blue','DisplayName','\omega_s');
hold on
plot(POT(:,9),POT(:,11), 'LineWidth', 2, 'Color', 'black','DisplayName','\omega_l');
plot(POT(:,9),POT(:,13), 'LineWidth', 2, 'Color', 'red','DisplayName','\omega_v');
plot(POT(:,9),2.*pi./POT(:,7), 'LineWidth', 2, 'Color', 'green','DisplayName','\omega');
xlabel('r1');
ylabel('\omega');
box('on');
grid on;
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)