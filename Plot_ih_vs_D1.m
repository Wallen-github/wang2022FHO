load('binast_d.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%i_h vs D_1\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_large(:,1)=table2array(binastd(table2array(binastd(:,4))>0.5,3));
data_large(:,2)=table2array(binastd(table2array(binastd(:,4))>0.5,25));
data_less(:,2)=table2array(binastd(table2array(binastd(:,4))<=0.5,25));
data_less(:,1)=table2array(binastd(table2array(binastd(:,4))<=0.5,3));
%scatter(table2array(binastd(:,3)),table2array(binastd(:,25)))
figure
scatter(data_large(:,1),data_large(:,2),'+','filled','blue')
hold on
scatter(data_less(:,1),data_less(:,2),'o','filled','black')
xlabel('D_A(km)');
ylabel('i_h (degree)');
box('on');
grid on;
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p_1 vs D_1\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_large(:,1)=table2array(binastd(table2array(binastd(:,4))>0.5,3));
data_large(:,2)=table2array(binastd(table2array(binastd(:,4))>0.5,8));
data_less(:,2)=table2array(binastd(table2array(binastd(:,4))<=0.5,8));
data_less(:,1)=table2array(binastd(table2array(binastd(:,4))<=0.5,3));
%scatter(table2array(binastd(:,3)),table2array(binastd(:,25)))
figure
scatter(data_large(:,1),data_large(:,2),'+','filled','blue')
hold on
scatter(data_less(:,1),data_less(:,2),'o','filled','black')
plot([min(table2array(binastd(:,3))) max(table2array(binastd(:,3)))],[2 2])
xlabel('D_A(km)');
ylabel('P_A (h)');
box('on');
grid on;
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p_1 vs D_1\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_large(:,1)=table2array(binastd(table2array(binastd(:,4))>0.5,3));
data_large(:,2)=table2array(binastd(table2array(binastd(:,4))>0.5,21));
data_less(:,2)=table2array(binastd(table2array(binastd(:,4))<=0.5,21));
data_less(:,1)=table2array(binastd(table2array(binastd(:,4))<=0.5,3));
%scatter(table2array(binastd(:,3)),table2array(binastd(:,25)))
figure
scatter(data_large(:,1),data_large(:,2),'+','filled','blue')
hold on
scatter(data_less(:,1),data_less(:,2),'o','filled','black')
xlabel('D_A(km)');
ylabel('\alpha_L');
box('on');
grid on;
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%p_2 vs ab_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_P2 = string(table2array(binastd(:,10)));
data_P2(~ismissing(data_P2)) = regexp(data_P2(~ismissing(data_P2)), '-?\d*\.?\d*', 'match');
data_P2 = double(string(data_P2));
data_Porb = string(table2array(binastd(:,9)));
data_Porb(data_Porb~='') = regexp(data_Porb(data_Porb~=''), '-?\d*\.?\d*', 'match');
data_Porb = double(string(data_Porb));
data_ab2 = string(table2array(binastd(:,16)));
data_ab2(~ismissing(data_ab2)) = regexp(data_ab2(~ismissing(data_ab2)), '-?\d*\.?\d*', 'match');
data_ab2 = double(string(data_ab2));

data_syn(:,1) = data_P2(abs(data_P2-data_Porb)<=5);
data_syn(:,2) = data_ab2(abs(data_P2-data_Porb)<=5);
data_asy(:,1) = data_P2(abs(data_P2-data_Porb)>5);
data_asy(:,2) = data_ab2(abs(data_P2-data_Porb)>5);
figure
scatter(data_P2,data_ab2,'filled','black')
hold on
scatter(data_syn(:,1),data_syn(:,2),'filled','black')
scatter(data_asy(:,1),data_asy(:,2),'filled','blue')
%plot([min(data_P2) max(data_P2)],[1.6 1.6])
%plot([20 20],[min(data_ab2) max(data_ab2)])
xlabel('P_B(h)');
ylabel('a/b_B');
box('on');
grid on;
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)