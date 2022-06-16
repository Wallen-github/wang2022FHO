DY = load('DeltaY.DAT');

figure
plot(DY(:,1),DY(:,2), 'LineWidth', 2, 'Color', 'blue','DisplayName','\omega_s');
hold on
plot(DY(:,1),DY(:,3), 'LineWidth', 2, 'Color', 'black','DisplayName','\omega_l');
xlabel('\delta Y');
ylabel('\delta Y');
box('on');
grid on;
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)