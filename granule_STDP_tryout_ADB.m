

%% plot both size distributions in one graph

figureHandle=figure;
set(gcf,'units','centimeters','position',[4,4,8,12.2])
load('kesten_workspace.mat')

%% Intrinsic dists
in_inds=[3,4,5];

all_lims=[-3,0.3 ; 0, 0.01 ; 0 ,1];


positions_X=([0.25, 2.25, 4.25, 6.25]);
positions_Y=([10.25, 8.25, 6.25, 4.25]);

labels={'Prop. stim.','Add. noise','Mult. noise'};

for ward=1:3
    ahMn = axes('Units','centimeters','Position',[positions_X(ward+1),positions_Y(ward),1.5,1.5]);
    pdSix = fitdist(M(:,in_inds(ward)),'Kernel');
    xgrid=linspace(all_lims(ward,1),all_lims(ward,2),100);
    ySix = pdf(pdSix,xgrid);
    area(xgrid,ySix,'EdgeColor','blue','LineWidth',1,'FaceColor','blue','FaceAlpha',0.25)
    set(gca,'LineWidth',0.75);
    set(gca,'Box','off')
    xlim([all_lims(ward,1),all_lims(ward,2)])
    set(gca,'YTick',[])
    ax=gca;
    ax.FontName='Helvetica';
    ax.FontSize=8;
    title(labels{ward})

    tprint('intrinsic properties','-jpg -tif -eps -png',[8,12.2])
   
end
pairs=[3,4;3,5;4,5];
inv_inds=[1,2;1,3;2,3];
idx=1;
colormap winter
for ii=1:2
    for jj=1:(ii)
      hMn = axes('Units','centimeters','Position',[positions_X(2+ii),positions_Y(jj),1.5,1.5]);  
      gridx1 = linspace(all_lims(inv_inds(idx,1),1),all_lims(inv_inds(idx,1),2),50);
      gridx2 = linspace(all_lims(inv_inds(idx,2),1),all_lims(inv_inds(idx,2),2),50);
      [x1,x2] = meshgrid(gridx1, gridx2);
      x1 = x1(:);
      x2 = x2(:);
      xi = [x1 x2];
      ksdensity([M(:,pairs(idx,1)),M(:,pairs(idx,2))],xi,'PlotFcn','surf')
      idx=idx+1;
      set(gca().Children,'EdgeColor','interp')
      view(hMn,[90 90]);
      grid(hMn,'off');
      ax=gca;
      ax.FontName='Helvetica';
      ax.FontSize=8;
      xlim([gridx1(1),gridx1(end)])
      ylim([gridx2(1),gridx2(end)])

      colormap(ax,'winter')
 
      set(gca,'XTick',[])
      set(gca,'YTick',[])
      tprint('intrinsic comparisons','-jpg -tif -eps -png',[8,12.2])

    end
end




%% Extrinsic dists
in_inds=[1,6,7];

all_lims=[100,500 ; 0, 5 ; 0 ,1];


positions_X=([0.25, 2.25, 4.25, 6.25]);
positions_Y=([9.95, 7.95, 5.95, 3.95]);

labels={'Tau_{STDP}','Pot. rate','Dep. rate'};

M(:,1)=10.^M(:,1); % Correct for exponentials
for ward=1:3
    ahMn = axes('Units','centimeters','Position',[positions_X(ward),positions_Y(ward),1.5,1.5]);
    pdSix = fitdist(M(:,in_inds(ward)),'Kernel');

    xgrid=linspace(all_lims(ward,1),all_lims(ward,2),100);
    ySix = pdf(pdSix,xgrid);
    area(xgrid,ySix,'EdgeColor','green','LineWidth',1,'FaceColor','green','FaceAlpha',0.25)
    set(gca,'LineWidth',0.75);
    set(gca,'Box','off')
    xlim([all_lims(ward,1),all_lims(ward,2)])
        set(gca,'YTick',[])
    ax=gca;
    ax.FontName='Helvetica';
    ax.FontSize=8;
    title(labels{ward})
    tprint('extrinsic properties','-jpg -tif -eps -png',[8,12.2])

   
end
pairs=[1,6;1,7;6,7];
inv_inds=[1,2;1,3;2,3];
idx=1;
for ii=1:2
    for jj=1:(ii)
      hMn = axes('Units','centimeters','Position',[positions_X(jj),positions_Y(ii+1),1.5,1.5]);  
      gridx1 = linspace(all_lims(inv_inds(idx,1),1),all_lims(inv_inds(idx,1),2),50);
      gridx2 = linspace(all_lims(inv_inds(idx,2),1),all_lims(inv_inds(idx,2),2),50);
      [x1,x2] = meshgrid(gridx1, gridx2);
      x1 = x1(:);
      x2 = x2(:);
      xi = [x1 x2];
      ksdensity([M(:,pairs(idx,1)),M(:,pairs(idx,2))],xi,'PlotFcn','surf')
      idx=idx+1;
      set(gca().Children,'EdgeColor','interp')
      view(hMn,[90 90]);
      grid(hMn,'off');
      ax=gca;
      ax.FontName='Helvetica';
      ax.FontSize=8;
      xlim([gridx1(1),gridx1(end)])
      ylim([gridx2(1),gridx2(end)])
      colormap(ax,'summer')
  
      set(gca,'XTick',[])
      set(gca,'YTick',[])
      tprint('extrinsic comparisons','-jpg -tif -eps -png',[8,12.2])

    end
end
%%
ahMn = axes('Units','centimeters','Position',[0.25,2.75,2.25,2.25]);
hold on
[x,c]=hist(D_silent(D_silent<1.5),30);
plot(c,x/trapz(c,x),'black')
iN_silent=N_silent(:);
iN_silent=iN_silent*mean(D_silent)/mean(iN_silent);
[X,C]=hist(iN_silent(iN_silent<1.5),30);
plot(C,X/trapz(C,X),'black','LineStyle','--')
set(gca,'LineWidth',0.75);
set(gca,'Box','off')
xlim([all_lims(ward,1),all_lims(ward,2)])
set(gca,'YTick',[])
ax=gca;
ax.FontName='Helvetica';
ax.FontSize=8;
title('Silent')
% tprint('Silent distribution','-jpg -tif -eps -png',[8,12.2])


ahMn = axes('Units','centimeters','Position',[2.75,2.75,2.25,2.25]);
hold on
plot(D_10(1,:),D_10(2,:),'black')
iN_10=N_10(:);
iN_10=iN_10*trapz(D_10(1,:),D_10(1,:).*D_10(2,:))/mean(iN_10);
[X,C]=hist(iN_10(iN_10<1.5),30);
plot(C,X/trapz(C,X),'black','LineStyle','--')
set(gca,'LineWidth',0.75);
set(gca,'Box','off')
xlim([all_lims(ward,1),all_lims(ward,2)])
set(gca,'YTick',[])
ax=gca;
ax.FontName='Helvetica';
ax.FontSize=8;
title('Control')
% tprint('control distribution','-jpg -tif -eps -png',[8,12.2])


ahMn = axes('Units','centimeters','Position',[5.5,2.75,2.25,2.25]);
hold on
plot(D_100(1,:),D_100(2,:),'black','DisplayName','Data')
iN_100=N_100(:);
iN_100=iN_100*trapz(D_100(1,:),D_100(1,:).*D_100(2,:))/mean(iN_100);
[X,C]=hist(iN_100(iN_100<1.5),30);
plot(C,X/trapz(C,X),'black','LineStyle','--','DisplayName','Model')
set(gca,'LineWidth',0.75);
set(gca,'Box','off')
xlim([all_lims(ward,1),all_lims(ward,2)])
set(gca,'YTick',[])
ax=gca;
ax.FontName='Helvetica';
ax.FontSize=8;
title('HFS')
leg=legend('boxoff');
leg.Units="centimeters";
leg.Position=[6.1,4.1,2,1];
leg.FontSize=8;
leg.FontName='Helvetica';
leg.ItemTokenSize=[20,8];
tprint('distributions','-jpg -tif -eps -png',[8,12.2])


tprint('Fig_8b','-SHR -jpg -tif -eps -png',[8,12.2])