%clear all;
%close all;
%FLI = load('RELATION.DAT');
%XY = load('RELATIONXY.DAT');
Numx = max(FLI(:,1));
Numy = max(FLI(:,2));
FLIMatrix = zeros(Numx,Numy);
TT = zeros(Numx,Numy);
Y = zeros(Numx,Numy);
X = zeros(Numx,Numy);
for i=1:Numx
    for j=1:Numy
        FLIMatrix(i,j)=FLI((i-1)*Numy+j,3);
        TT(i,j)=FLI((i-1)*Numy+j,4);
        Y(i,j)=XY((i-1)*Numy+j,3);
        X(i,j)=XY((i-1)*Numy+j,4);
    end
end

% for i=1:length(FLI)
%     FLIMatrix(int8(FLI(i,1)),int8(FLI(i,2)))=FLI(i,3);
%     TT(int8(FLI(i,1)),int8(FLI(i,2)))=FLI(i,4);
%     Y(int8(XY(i,1)),int8(XY(i,2)))=XY(i,3);
%     X(int8(XY(i,1)),int8(XY(i,2)))=XY(i,4);
% end
%FLIMatrix(FLIMatrix==0) = NaN;
% FLIMatrix(FLIMatrix>299) = NaN;

figure;
contourf(X,Y,FLIMatrix,10)
%xticks([0 Numx/6 Numx*2/6 Numx*3/6 Numx*4/6 Numx*5/6 Numx])
%xticklabels({X(1,1),X(1,int8(Numx/6)),X(1,int8(Numx*2/6)),X(1,int8(Numx*3/6)),X(1,int8(Numx*4/6)),X(1,int8(Numx*5/6)),X(1,int8(Numx))})
%yticks([0 Numy/6 Numy*2/6 Numy*3/6 Numy*4/6 Numy*5/6 Numy])
%yticklabels({Y(1,1),Y(int8(Numy/6),1),Y(int8(Numy*2/6),1),Y(int8(Numy*3/6),1),Y(int8(Numy*4/6),1),Y(int8(Numy*5/6),1),Y(int8(Numy),1)})
xlabel('$r^\prime$ (unit length)','interpreter','latex');
ylabel('$r_\odot$ (AU)','interpreter','latex');
shading flat
colorbar
set(gca,'Clim',[0 40])
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

figure;
surf(X,Y,FLIMatrix);
shading flat
colorbar
xlabel('$r^\prime$ (unit length)','interpreter','latex');
ylabel('$r_\odot$ (AU)','interpreter','latex');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)
%set(gca,'Xlim',[5 14])
%set(gca,'Zlim',[25 35])
set(gca,'Clim',[0 40])


