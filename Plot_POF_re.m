load('POF.mat')
POT = POT1_re;


figure
subplot(2,1,1)
plot(POT(:,8),POT(:,2), 'LineWidth', 2, 'Color', 'black');
ylabel('y');
box('on');
grid on;
set(gca,'xticklabel',[]);
set(gca,'xlim',[15 160]);
%set(gca,'ylim',[0.97 1.02]);
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
set(gca,'xlim',[15 160]);
set(gca,'ylim',[0.5 1.5]);
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2,'position',[0.15 0.15 0.75 0.4])