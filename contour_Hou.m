%%
clear all;
clc;
%% 读取数据
idx1= 16;
idx2= 12;
data1=importdata('contour.txt');
xaxislim=0.15;
yaxislim=0.5;
xchar='$TOF$';
ychar='$\Delta T$';
filename='figure1.fig';

%% 写入格点
M=301;
N=31;

DRIFT=zeros(M,N);

for j=1:M
    for i=1:N
        
    if(data1((j-1)*N+i,3)==0)
        data1((j-1)*N+i,3)=NaN;
    end
%%    DRIFT(j,i)=log10(data1((j-1)*N+i,3));
       DRIFT(j,i)=data1((j-1)*N+i,3);
    end
end

x=1:N;
y=1:M;
[X,Y]=meshgrid(x,y);

%% 作图
figure;
surf(X,Y,DRIFT);
shading flat; 
% shading interp; 
colorbar;
axis([0 30 0 300]);
box on;
colormap JET;
view([0 0 1]); 

title('$coutour map$','interpreter','latex');
ylabel(ychar,'interpreter','latex');
xlabel(xchar,'interpreter','latex');
set(gca,'xtick',[0,10,20,30],'xticklabel',0:1:3);
set(gca,'ytick',[0,100,200,300],'yticklabel',0:10:30);
