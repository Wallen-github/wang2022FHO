%close all;
%clear all;
POT = load('POF.DAT');
%load('POF.mat')
%POT = POT1_re;

figure
plot(POT(:,7),POT(:,2), 'LineWidth', 2, 'Color', 'black');
xlabel('Tp');
ylabel('S');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

figure
plot(POT(:,9),POT(:,2), 'LineWidth', 2, 'Color', 'black');
xlabel('axis(m)');
ylabel('S');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)


figure
subplot(2,1,1)
for i=1:length(POT(:,1))
    if POT(i,10) == 1
        scatter(POT(i,8),POT(i,2), '.','blue');
        hold on
        scatter(POT(i,8),POT(i,14), '.','green');
        scatter(POT(i,8),POT(i,15), '.','green');
    elseif POT(i,10) == 0
        scatter(POT(i,8),POT(i,2), '.','magenta');
        hold on
        scatter(POT(i,8),POT(i,14), '.','green');
        scatter(POT(i,8),POT(i,15), '.','green');
    end
end
%plot([min(POT2(i,2)) max(POT2(i,2))],[0.909090909090909 0.909090909090909],'LineWidth', 2,'Color','red','LineStyle','--','DisplayName','b_A');
xlabel('Tn');
ylabel('Yp');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

subplot(2,1,2)
plot(POT(:,8),POT(:,11), 'LineWidth', 2, 'Color', 'black','DisplayName','\omega_l');
hold on
plot(POT(:,8),POT(:,12), 'LineWidth', 2, 'Color', 'blue','DisplayName','\omega_s');
%plot(POT(:,8),POT(:,13), 'LineWidth', 2, 'Color', 'red','DisplayName','\omega_v');
plot(POT(:,8),2.*pi./POT(:,7), 'LineWidth', 2, 'Color', 'green','DisplayName','\omega_p');
xlabel('Tn');
ylabel('\omega');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)
