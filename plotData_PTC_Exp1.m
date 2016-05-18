%%% Exp. 1 in PTC article
%% off-site listening? 
% plot CURVES WITHOUT bias
% Aug.14
clear all; clc;
figure(152)
dat = '26-Mar-2015';
% Subj.3
clearvars -except HW dat; load(strcat('MedEl_S3_SM_M2_',dat));
subplot(2,4,1);
%     plot([-10 5],[ThM2 ThM2],'--','color',[.5 .5 .5]);hold on;
    plot( 0:4,Tm2,'vk',0:4,Cm2,'^k','markersize',5,'markerfacecolor','k'); hold on;
    plot(-4:0,Tm2(end:-1:1),'vk',-4:0,Cm2(end:-1:1),'^k','markersize',5,'markerfacecolor','k');
    plot( 0:4,Tm1(7:11),'>m',-4:0,Tm1(3:1:7),'<g','markersize',5,'markerfacecolor','b','color','b'); % plot Thresholds for individual maskers
%     plot([-4],Cm1([3]),'*g',[4],Cm1([11]),'*m','color',grey,'MarkerSize',8,'linewidth',1);   
    errorbar(0:2,M1(4:6),e1(4:6),'o-m','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-2:0,M1(2:4),e1(2:4),'o-g','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-4:4,[M2(end:-1:2) M2],[e2(end:-1:2) e2],'o-k','linewidth',1.5,'markerfacecolor','k'); 
    hold off;
    ymsi=.75; x=0:0.1:3; 
    apical  = fliplr((M1(1:4)-min(M1(1:4)))/(max(M1(1:4))-min(M1(1:4)))); 
    yintB=interp1(min(x):max(x),apical,x); HWb=x(argmin(abs(yintB-ymsi)));
    basal = ((M1(4:7)-min(M1(4:7)))/(max(M1(4:7))-min(M1(4:7)))); 
    yintA=interp1(min(x):max(x),basal,x); HWa=x(argmin(abs(yintA-ymsi)));
    dual   = ((M2(1:4)-min(M2(1:4)))/(max(M2(1:4))-min(M2(1:4)))); 
    yintD=interp1(min(x):max(x),dual,x); HWd=x(argmin(abs(yintD-ymsi)));
    HW(1,:)=[HWb HWa HWd];
    
    ylabel({'MLTs for dual and single masker [dB]'}); xlabel('\Deltax');title('ME1'); 
    axis([-5 5 47 62]);      set(gca,'XTick',-4:4,'YTick',48:3:62); 
      
% Subj.4
clearvars -except HW dat;  load(strcat('MedEl_S4_AL_M2_',dat));
subplot(2,4,2);
%     plot([-10 5],[ThM2 ThM2],'--','color',[.5 .5 .5]);hold on;
    plot(0:4,Tm2,'vk',0:4,Cm2,'^k','markersize',5,'markerfacecolor','k');hold on;
    plot(-4:0,Tm2(end:-1:1),'vk',-4:0,Cm2(end:-1:1),'^k','markersize',5,'markerfacecolor','k');
    plot(0:4,Tm1(7:11),'>m',-4:0,Tm1(3:7),'<g','markersize',5,'markerfacecolor','b','color','b'); 
    errorbar(0:4,M1(5:end),e1(5:end),'o-m','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-4:0,M1(1:5),e1(1:5),'o-g','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-4:4,[M2(end:-1:2) M2],[e2(end:-1:2) e2],'o-k','linewidth',1.5,'markerfacecolor','k'); 
    ymsi=.75; x=0:0.1:4; 
    apical  = fliplr((M1(1:5)-min(M1(1:5)))/(max(M1(1:5))-min(M1(1:5)))); 
    yintB=interp1(min(x):max(x),apical,x); HWb=x(argmin(abs(yintB-ymsi)));
    basal = ((M1(5:9)-min(M1(5:9)))/(max(M1(5:9))-min(M1(5:9)))); 
    yintA=interp1(min(x):max(x),basal,x); HWa=x(argmin(abs(yintA-ymsi)));
    dual   = ((M2(1:5)-min(M2(1:5)))/(max(M2(1:5))-min(M2(1:5)))); 
    yintD=interp1(min(x):max(x),dual,x); HWd=x(argmin(abs(yintD-ymsi)));
    HW(2,:)=[HWb HWa HWd];
    hold off; 
    set(gca,'XTick',-4:4,'YTick',48:3:62); xlabel('\Deltax'); title('ME2'); 
    axis([-5 5 47 62]);
     
   
    
% Subj.8
clearvars -except HW dat;  load(strcat('MedEl_S8_AM_M2_',dat));
subplot(2,4,3);%%% simultaneous dual
%     plot([-10 5],[ThM2 ThM2],'--','color',[.5 .5 .5]);hold on;  
    plot(0:4,Tm2,'vk',0:4,Cm2,'^k','markersize',5,'markerfacecolor','k');hold on;
    plot(-4:0,Tm2(end:-1:1),'vk',-4:0,Cm2(end:-1:1),'^k','markersize',5,'markerfacecolor','k');
    plot(0:4,Tm1(7:11),'>m',-4:0,Tm1(3:7),'<g','markersize',5,'markerfacecolor','b','color','b');  
    plot(4,Cm2(end),'*','color','k','MarkerSize',8,'linewidth',1);
    errorbar(0:4,M1(5:end),e1(5:end),'o-m','linewidth',1.5,'color','b','markerfacecolor','w'); 
    errorbar(-4:0,M1(1:5),e1(1:5),'o-g','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-4:4,[M2(end:-1:2) M2],[e2(end:-1:2) e2],'o-k','linewidth',1.5,'markerfacecolor','k'); 
    ymsi=.75; x=0:0.1:4; 
    apical  = fliplr((M1(1:5)-min(M1(1:5)))/(max(M1(1:5))-min(M1(1:5)))); 
    yintB=interp1(min(x):max(x),apical,x); HWb=x(argmin(abs(yintB-ymsi)));
    basal = ((M1(5:9)-min(M1(5:9)))/(max(M1(5:9))-min(M1(5:9)))); 
    yintA=interp1(min(x):max(x),basal,x); HWa=x(argmin(abs(yintA-ymsi)));
    dual   = ((M2(1:5)-min(M2(1:5)))/(max(M2(1:5))-min(M2(1:5)))); 
    yintD=interp1(min(x):max(x),dual,x); HWd=x(argmin(abs(yintD-ymsi)));
    HW(3,:)=[HWb HWa HWd];
    hold off; 
    set(gca,'XTick',-4:4,'YTick',48:3:62); xlabel('\Deltax'); title('ME4'); 
    axis([-5 5 47 62]);

% Subj.9
clearvars -except HW dat;  load(strcat('MedEl_S9_CL_M2_',dat));
subplot(2,4,4);
%     plot([-10 5],[ThM2 ThM2],'--','color',[.5 .5 .5]);
    plot(0:4,Tm2,'vk',0:4,Cm2,'^k','markersize',5,'markerfacecolor','k');hold on;
    plot(-4:0,Tm2(end:-1:1),'vk',-4:0,Cm2(end:-1:1),'^k','markersize',5,'markerfacecolor','k');
    plot(0:4,Tm1((7:11)-2),'>m',-4:0,Tm1((3:7)-2),'<g','markersize',5,'markerfacecolor','b','color','b'); 
    plot(4,Cm2(end),'*','color','k','MarkerSize',8,'linewidth',1.5); 
    errorbar(0:2,M1(3:5),e1(3:5),'o-m','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-2:0,M1(1:3),e1(1:3),'o-g','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-4:4,[M2(end:-1:2) M2],[e2(end:-1:2) e2],'o-k','linewidth',1.5,'markerfacecolor','k'); 
    ymsi=.75; x=0:0.1:2; 
    apical  = fliplr((M1(1:3)-min(M1(1:3)))/(max(M1(1:3))-min(M1(1:3)))); 
    yintB=interp1(min(x):max(x),apical,x); HWb=x(argmin(abs(yintB-ymsi)));
    basal = ((M1(3:5)-min(M1(3:5)))/(max(M1(3:5))-min(M1(3:5)))); 
    yintA=interp1(min(x):max(x),basal,x); HWa=x(argmin(abs(yintA-ymsi)));
    dual   = ((M2(1:3)-min(M2(1:3)))/(max(M2(1:3))-min(M2(1:3)))); 
    yintD=interp1(min(x):max(x),dual,x); HWd=x(argmin(abs(yintD-ymsi)));
    HW(4,:)=[HWb HWa HWd];
    hold off; 
    set(gca,'XTick',-4:4,'YTick',48:3:62); xlabel('\Deltax'); title('ME5'); 
    axis([-5 5 47 62]);

%%% AB 
% Subj.1
clearvars -except HW dat;  load(strcat('AB_S1_AM_M2_',dat));
subplot(2,4,5);    
%     plot([-1 5],[ThM2 ThM2],'--','color',[.5 .5 .5]);   hold on
    plot(0:4,Tm2,'vk',0:4,Cm2,'^k','markersize',5,'markerfacecolor','k'); hold on;
    plot(-4:0,Tm2(end:-1:1),'vk',-4:0,Cm2(end:-1:1),'^k','markersize',5,'markerfacecolor','k');
    plot(0:2,Tm1(3:end),'>m',-2:0,Tm1([1 2 3]),'<g','markersize',5,'markerfacecolor','b','color','b');
    errorbar(0:2,M1(3:end),e1(3:end),'o-m','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-2:0, M1([1 2 3]),e1(1:3),'o-g','linewidth',1.5,'color','b','markerfacecolor','w'); 
    errorbar(0:4,M2,e2,'o-k','linewidth',2,'markerfacecolor','k');
    errorbar(-4:0,M2(end:-1:1),e2(end:-1:1),'o-k','linewidth',2,'markerfacecolor','k');    
    ymsi=.75; x=0:0.1:2; 
    apical  = fliplr((M1(1:3)-min(M1(1:3)))/(max(M1(1:3))-min(M1(1:3)))); 
    yintB=interp1(min(x):max(x),apical,x); HWb=x(argmin(abs(yintB-ymsi)));
    basal = ((M1(3:5)-min(M1(3:5)))/(max(M1(3:5))-min(M1(3:5)))); 
    yintA=interp1(min(x):max(x),basal,x); HWa=x(argmin(abs(yintA-ymsi)));
    dual   = ((M2(1:3)-min(M2(1:3)))/(max(M2(1:3))-min(M2(1:3)))); 
    yintD=interp1(min(x):max(x),dual,x); HWd=x(argmin(abs(yintD-ymsi)));
    HW(5,:)=[HWb HWa HWd];
    hold off;
    set(gca,'XTick',-4:4,'YTick',45:3:57); ylabel({'MLTs for dual and single masker [dB]'}); xlabel('\Deltax'); title('AB1'); 
    axis([-5 5 44 58]);     
    

% Subj.2
clearvars -except HW dat;  load(strcat('AB_S2_DM_M2_',dat));
subplot(2,4,6);
%     plot([-1 5],[ThM2 ThM2],'--','color',[.5 .5 .5]);
    plot(0:4,Tm2,'vk',0:4,Cm2,'^k','markersize',5,'markerfacecolor','k');hold on;
    plot(-4:0,Tm2(end:-1:1),'vk',-4:0,Cm2(end:-1:1),'^k','markersize',5,'markerfacecolor','k');
    plot(0:3,Tm1(4:end),'>m',-3:0,Tm1(1:4),'<g','markersize',5,'markerfacecolor','b','color','b');
    errorbar(0:3,M1(4:end),e1(4:end),'o-m','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(-3:0,M1(1:4),e1(1:4),'o-g','linewidth',1.5,'color','b','markerfacecolor','w');
    errorbar(0:4,M2,e2,'o-k','linewidth',1.5,'markerfacecolor','k'); 
    errorbar(-4:0,M2(end:-1:1),e2(end:-1:1),'o-k','linewidth',2,'markerfacecolor','k');
    ymsi=.75; x=0:0.1:3; 
    apical  = fliplr((M1(1:4)-min(M1(1:4)))/(max(M1(1:4))-min(M1(1:4)))); 
    yintB=interp1(min(x):max(x),apical,x); HWb=x(argmin(abs(yintB-ymsi)));
    basal = ((M1(4:7)-min(M1(4:7)))/(max(M1(4:7))-min(M1(4:7)))); 
    yintA=interp1(min(x):max(x),basal,x); HWa=x(argmin(abs(yintA-ymsi)));
    dual   = ((M2(1:4)-min(M2(1:4)))/(max(M2(1:4))-min(M2(1:4)))); 
    yintD=interp1(min(x):max(x),dual,x); HWd=x(argmin(abs(yintD-ymsi)));
    HW(6,:)=[HWb HWa HWd];
    hold off; 
    set(gca,'XTick',-4:4,'YTick',45:3:57); xlabel('\Deltax'); title('AB2'); 
    axis([-5 5 44 58]);
   
subplot(2,4,7)
plot(0,0,'ok','markerfacecolor','k','linewidth',1.5); hold on;
plot(0,0,'o', 'markerfacecolor','w','color','b','linewidth',1.5);
plot(0,0,'kv',0,0,'k^','markerfacecolor','k','linewidth',1.5,'markersize',5);
plot(0,0,'b>',0,0,'b<','markerfacecolor','b','linewidth',1.5,'markersize',5); 
plot(0,0,'--','color',[.5 .5 .5]); hold off; xlim([3 4]); axis off;
legend('MLT dual masker', 'MLT single','T dual','C dual','T single basal','T single apical','location','North');

set(gcf,'position',[70 10 1100 520]);
% print(gcf, '-dpng', '-r200', 'Fig_4_PTC.png');

%% measure stats on widths
[h, p] = ttest(HW(:,2),HW(:,1)); res0(1) = p; % apical Vs basal
[h, p] = ttest(HW(:,2),HW(:,3)); res0(2) = p; % apical Vs dual
[h, p] = ttest(HW(:,1),HW(:,3)); res0(3) = p; % basal Vs dual
[h, p] = ttest(nanmin(HW(:,[1 2])')',HW(:,3)); res0(4) = p; % least effective Vs dual
[h, p] = ttest(nanmax(HW(:,[1 2])')',HW(:,3)); res0(5) = p; % most effective Vs dual
[h, p, ci, stats] = ttest(sum(HW(:,[1 2])')', 2*HW(:,3)); res0(7) = p; % sum single Vs dual (Dingemanse's method)


%% off-site listening? plot CURVES
% %%% MedElAug.14
clear all; clc; closefigs;
flag_do_stats = false;
figure(149)
D=zeros(5,6); S=zeros(5,6); 

% Subj.3
clearvars -except D S p_slope P flag_do_stats; dat = '26-Mar-2015';  load(strcat('MedEl_S3_SM_M2_',dat));
subplot(2,4,1);%%% simultaneous dual and single
    IDme = [4 3 2]; IDle = [4 5 6]; dx=0:2;
    res_M1mem = 20*log10( 10.^(b1(IDme,:)./20)) - (M1(4)- M2(1));%- (10.^(M1(4)/20) - 10.^(M2(1)/20)) );
    res_M1lem = 20*log10( 10.^(b1(IDle,:)./20)) - (M1(4)- M2(1));%- (10.^(M1(4)/20) - 10.^(M2(1)/20)) );
    M1mem     = nanmean(res_M1mem,2); M1memSD     = nanstd(res_M1mem,[],2);
    M1lem     = nanmean(res_M1lem,2); M1lemSD     = nanstd(res_M1lem,[],2);
    plot(1,M1lem(2),'<',2,M1lem(3),'>','linewidth',2,'color',grey,'markerfacecolor',grey,'markersize',5);hold on;
    errorbar(dx-.05,M2(1:3),e2(1:3),'o-k','linewidth',1.5,'markerfacecolor','k'); 
    errorbar(dx+.05,M1mem,  M1memSD,'o-r','linewidth',1.5,'markerfacecolor','w'); 
    D(1:3,1)=M2(1:3); S(1:3,1)=M1mem; 
    for n=1:length(IDme), [tmp1, tmp2] = ttest2(res_M1mem(n,:),b2(n,:)); 
        s_vals(n) = tmp1; p_vals(n) = tmp2; end;
    [p2, k2] = polyfit(dx(1:3),M2(1:3),1); [p1, k1] = polyfit(dx(1:3),M1mem(1:3)',1);
    Z=(p1(1)-p2(1)) / sqrt(k1.normr^2 + k2.normr^2);
    p_slope(1) = normcdf(-abs(Z),0,1);
    plot(dx(logical(s_vals)),61.5*ones(size(dx(logical(s_vals)))),'*k')
    hold off;
    set(gca,'XTick',-0:4); title('ME1'); 
    axis([-1 5 50 60]);     
    ylabel({'MLTs for dual and single masker [dB]'}); xlabel('\Deltax');
%     y  = [clm((b2(1:3,:))')'  clm(res_M1mem')']';
%     g1 = [1*ones(1,4*3)             2*ones(1,4*3)]; 
%     g2 = [clm(((0:2)'*ones(1,4))')' clm(((0:2)'*ones(1,4))')']; 
    y  = [clm((b2(2:3,:))')'  clm(res_M1mem(2:end,:)')']'; %MLTs
    g1 = [1*ones(1,4*2)             2*ones(1,4*2)];              %type
    g2 = [clm(((1:2)'*ones(1,4))')' clm(((1:2)'*ones(1,4))')'];  %dx
    group={g1 g2}; if flag_do_stats, [p,table,stats,terms] = anovan(y,group,'model',[1 0;0 1;1 1]); end;
    P(1,:)=p;

% Subj.4
clearvars -except D S p_slope P flag_do_stats; dat = '26-Mar-2015';  load(strcat('MedEl_S4_AL_M2_',dat));
subplot(2,4,2);%%% simultaneous dual
    IDme = [5 4 3 2 1]; IDle = [5 6 7 8 9];  dx=0:4;
    res_M1mem = 20*log10( 10.^(b1(IDme,:)./20)) - (M1(5)- M2(1));% - (10.^(M1(4)/20) - 10.^(M2(1)/20)) );
    res_M1lem = 20*log10( 10.^(b1(IDle,:)./20)) - (M1(5)- M2(1));% - (10.^(M1(4)/20) - 10.^(M2(1)/20)) );
    M1mem     = nanmean(res_M1mem,2); M1memSD     = nanstd(res_M1mem,[],2);
    M1lem     = nanmean(res_M1lem,2); M1lemSD     = nanstd(res_M1lem,[],2);
    plot(1,M1lem(2),'<',2,M1lem(3),'<',3,M1lem(4),'>','linewidth',2,'color',grey,'markerfacecolor',grey,'markersize',5);hold on;
    errorbar(dx-.05,M2(1:5),e2(1:5),'o-k','linewidth',1.5,'markerfacecolor','k'); 
    errorbar(dx+.05,M1mem,  M1memSD,'o-r','linewidth',1.5,'markerfacecolor','w'); 
%     D(1:4,2)=M2(1:4); S(1:4,2)=min([M1(4:7)-(M1(4)-M2(1)); M1(4:-1:1)-(M1(4)-M2(1))]); 
    D(:,2)=M2; S(:,2)=min([M1(5:9)-(M1(5)-M2(1)); M1(5:-1:1)-(M1(5)-M2(1))]);     
    for n=1:length(IDme), [tmp1, tmp2] = ttest2(res_M1mem(n,:),b2(n,:)); 
        s_vals(n) = tmp1; p_vals(n) = tmp2; end;
    plot(dx(logical(s_vals)),61.5*ones(size(dx(logical(s_vals)))),'*k')
    [p2, k2] = polyfit(dx(1:end),M2(1:length(IDme)),1); [p1, k1] = polyfit(dx(1:end),M1mem(1:end)',1);
    Z=(p1(1)-p2(1)) / sqrt(k1.normr^2 + k2.normr^2);
    p_slope(2) = normcdf(-abs(Z),0,1);
    hold off; 
    set(gca,'XTick',-0:4); title('ME2'); 
    axis([-1 5 50 60]);xlabel('\Deltax');
%     y  = [clm((b2-M2(1))')'  clm(res_M1mem')']';
%     g1 = [1*ones(1,4*5)             2*ones(1,3*5)]; 
%     g2 = [clm(((0:4)'*ones(1,4))')' clm(((0:4)'*ones(1,3))')']; 
    y  = [clm((b2(2:end,:))')'      clm(res_M1mem(2:end,:)')']';
    g1 = [1*ones(1,4*4)             2*ones(1,4*4)]; 
    g2 = [clm(((1:4)'*ones(1,4))')' clm(((1:4)'*ones(1,4))')']; 
    group={g1 g2}; if flag_do_stats, [p,table,stats,terms] = anovan(y,group,'model',[1 0;0 1;1 1]); end;
    P(2,:)=p;

% Subj.8
clearvars -except D S p_slope P flag_do_stats; dat = '26-Mar-2015';  load(strcat('MedEl_S8_AM_M2_',dat));
subplot(2,4,3);%%% simultaneous dual
    IDme = [5 6 3 2 1]; IDle = [5 4 7 8 1];  dx=0:4;
    res_M1mem = 20*log10( 10.^(b1(IDme,:)./20)) - (M1(5)- M2(1));% - (10.^(M1(5)/20) - 10.^(M2(1)/20)) );
    res_M1lem = 20*log10( 10.^(b1(IDle,:)./20)) - (M1(5)- M2(1));% - (10.^(M1(5)/20) - 10.^(M2(1)/20)) );
    M1mem     = nanmean(res_M1mem,2); M1memSD     = nanstd(res_M1mem,[],2);
    M1lem     = nanmean(res_M1lem,2); M1lemSD     = nanstd(res_M1lem,[],2);
    plot(1,M1lem(2),'>',2,M1lem(3),'<',3,M1lem(4),'<',4,M1lem(5),'<','linewidth',2,'color',grey,'markerfacecolor',grey,'markersize',5);hold on;
    errorbar(dx-0.05,M2,e2,'o-k','linewidth',1.5,'markerfacecolor','k'); 
    errorbar(dx+0.05,M1mem,  M1memSD,'o-r','linewidth',1.5,'markerfacecolor','w'); 
    D(:,3)=M2; S(:,3)=min([M1(5:9)-(M1(5)-M2(1)); M1(5:-1:1)-(M1(5)-M2(1))]); 
    for n=1:length(IDme), [tmp1, tmp2] = ttest2(res_M1mem(n,:),b2(n,:)); 
        s_vals(n) = tmp1; p_vals(n) = tmp2; end;
    plot(dx(logical(s_vals)),61.5*ones(size(dx(logical(s_vals)))),'*k')
    [p2, k2] = polyfit(dx,M2,1); [p1, k1] = polyfit(dx,M1mem',1);
    Z=(p1(1)-p2(1)) / sqrt(k1.normr^2 + k2.normr^2);
    p_slope(3) = normcdf(-abs(Z),0,1);
    hold off; 
    set(gca,'XTick',-0:4); title('ME4'); 
    axis([-1 5 50 60]);xlabel('\Deltax');
%     y  = [clm((b2-M2(1))')'  clm(res_M1mem')']';
%     g1 = [1*ones(1,4*5)             2*ones(1,2*5)]; 
%     g2 = [clm(((0:4)'*ones(1,4))')' clm(((0:4)'*ones(1,2))')']; 
    y  = [clm((b2(2:end,:))')'      clm(res_M1mem(2:end,:)')']';
    g1 = [1*ones(1,4*4)             2*ones(1,2*4)]; 
    g2 = [clm(((1:4)'*ones(1,4))')' clm(((1:4)'*ones(1,2))')']; 
    group={g1 g2}; if flag_do_stats, [p,table,stats,terms] = anovan(y,group,'model',[1 0;0 1;1 1]); end;
    P(3,:)=p;
% Subj.9
clearvars -except D S p_slope P flag_do_stats; dat = '26-Mar-2015';  load(strcat('MedEl_S9_CL_M2_',dat));
subplot(2,4,4);%%% simultaneous dual
    IDme = [3 4 5]; IDle = [3 2 1];  dx=0:2;
    res_M1mem = 20*log10( 10.^(b1(IDme,:)./20)) - (M1(3)- M2(1));
    res_M1lem = 20*log10( 10.^(b1(IDle,:)./20)) - (M1(3)- M2(1));
    M1mem     = nanmean(res_M1mem,2); M1memSD     = nanstd(res_M1mem,[],2);
    M1lem     = nanmean(res_M1lem,2); M1lemSD     = nanstd(res_M1lem,[],2);
    plot(1,M1lem(2),'<',2,M1lem(3),'<','linewidth',2,'color',grey,'markerfacecolor',grey,'markersize',5);hold on;
    errorbar(dx-.05,M2(1:3),e2(1:3),'o-k','linewidth',1.5,'markerfacecolor','k'); 
    errorbar(dx+.05,M1mem,  M1memSD,'o-r','linewidth',1.5,'markerfacecolor','w'); 
    D(1:3,4)=M2(1:3); S(1:3,4)=min([M1(3:5)-(M1(3)-M2(1)); M1(3:-1:1)-(M1(3)-M2(1))]); 
    for n=1:length(IDme), [tmp1, tmp2] = ttest2(res_M1mem(n,:),b2(n,:)); 
        s_vals(n) = tmp1; p_vals(n) = tmp2; end;
    plot(dx(logical(s_vals)),61.5*ones(size(dx(logical(s_vals)))),'*k')
    [p2, k2] = polyfit(dx,M2(1:length(IDme)),1); [p1, k1] = polyfit(dx,M1mem(1:end)',1);
    Z=(p1(1)-p2(1)) / sqrt(k1.normr^2 + k2.normr^2);
    p_slope(4) = normcdf(-abs(Z),0,1);
    hold off; 
    set(gca,'XTick',-0:4); title('ME5'); 
    axis([-1 5 50 60]);xlabel('\Deltax');
%     y  = [clm((b2(1:3,:))')'  clm(res_M1mem')']';
%     g1 = [1*ones(1,4*3)             2*ones(1,4*3)]; 
%     g2 = [clm(((0:2)'*ones(1,4))')' clm(((0:2)'*ones(1,4))')']; 
    y  = [clm((b2(2:3,:))')'        clm(res_M1mem(2:end,:)')']';
    g1 = [1*ones(1,4*2)             2*ones(1,4*2)]; 
    g2 = [clm(((1:2)'*ones(1,4))')' clm(((1:2)'*ones(1,4))')']; 
    group={g1 g2}; if flag_do_stats, [p,table,stats,terms] = anovan(y,group,'model',[1 0;0 1;1 1]); end;
    P(4,:)=p;
%%% AB 
% Subj.1
clearvars -except D S p_slope P flag_do_stats; dat = '26-Mar-2015';  load(strcat('AB_S1_AM_M2_',dat));
subplot(2,4,5);%%% simultaneous dual and single
    IDme = [3 2 1]; IDle = [3 4 1];  dx=0:2;
    res_M1mem = 20*log10( 10.^(b1(IDme,:)./20)) - (M1(3)- M2(1));% - (10.^(M1(3)/20) - 10.^(M2(1)/20)) );
    res_M1lem = 20*log10( 10.^(b1(IDle,:)./20)) - (M1(3)- M2(1));% - (10.^(M1(3)/20) - 10.^(M2(1)/20)) );
    M1mem     = nanmean(res_M1mem,2); M1memSD     = nanstd(res_M1mem,[],2);
    M1lem     = nanmean(res_M1lem,2); M1lemSD     = nanstd(res_M1lem,[],2);
    plot(1,M1lem(2),'<',2,M1lem(3),'<','linewidth',2,'color',grey,'markerfacecolor',grey,'markersize',5);hold on;
    errorbar(dx-.05,M2(1:3),e2(1:3),'o-k','linewidth',1.5,'markerfacecolor','k'); 
    errorbar(dx+.05,M1mem,  M1memSD,'o-r','linewidth',1.5,'markerfacecolor','w'); 
    D(1:3,5)=M2(1:3); S(1:3,5)=min([M1(3:5)-(M1(3)-M2(1)); M1(3:-1:1)-(M1(3)-M2(1))]); 
    for n=1:length(IDme), [tmp1, tmp2] = ttest2(res_M1mem(n,:),b2(n,:)); 
        s_vals(n) = tmp1; p_vals(n) = tmp2; end;
%     plot(dx(logical(s_vals)),55*ones(size(dx(logical(s_vals)))),'*k')
    [p2, k2] = polyfit(dx,M2(1:3),1); [p1, k1] = polyfit(dx,M1mem',1);
    Z=(p1(1)-p2(1)) / sqrt(k1.normr^2 + k2.normr^2);
    p_slope(5) = normcdf(-abs(Z),0,1);
    hold off;
    set(gca,'XTick',0:4,'Ytick',45:2:55); title('AB1'); 
    axis([-1 5 45 55]);     
    ylabel({'MLTs for dual and single masker [dB]'}); xlabel('\Deltax');
    y  = [clm((b2(2:3,:))')'        clm(res_M1mem(2:3,:)')']';
    g1 = [1*ones(1,4*1)             2*ones(1,4*2)]; %type
    g2 = [clm(((1:2)'*ones(1,2))')' clm(((1:2)'*ones(1,4))')']; %dx
    group={g1 g2}; if flag_do_stats, [p,table,stats,terms] = anovan(y,group,'model',[1 0;0 1;1 1]); end;
    P(5,:)=p;
% Subj.2
clearvars -except D S p_slope P flag_do_stats; dat = '26-Mar-2015';  load(strcat('AB_S2_DM_M2_',dat));
subplot(2,4,6);%%% simultaneous dual
    IDme = [4 3 2 7]; IDle = [4 5 6 1];  dx=0:3;
    res_M1mem = 20*log10( 10.^(b1(IDme,:)./20)) - (M1(4)- M2(1));% - (10.^(M1(4)/20) - 10.^(M2(1)/20)) );
    res_M1lem = 20*log10( 10.^(b1(IDle,:)./20)) - (M1(4)- M2(1));% - (10.^(M1(4)/20) - 10.^(M2(1)/20)) );
    M1mem     = nanmean(res_M1mem,2); M1memSD     = nanstd(res_M1mem,[],2);
    M1lem     = nanmean(res_M1lem,2); M1lemSD     = nanstd(res_M1lem,[],2);
    plot(1,M1lem(2),'<',2,M1lem(3),'<',3,M1lem(4),'>','linewidth',2,'color',grey,'markerfacecolor',grey,'markersize',5);hold on;
    errorbar(dx-.05,M2(1:4),e2(1:4),'o-k','linewidth',1.5,'markerfacecolor','k'); 
    errorbar(dx+.05,M1mem,  M1memSD,'o-r','linewidth',1.5,'markerfacecolor','w'); 
    D(1:4,6)=M2(1:4); S(1:4,6)=min([M1(4:7)-(M1(4)-M2(1)); M1(4:-1:1)-(M1(4)-M2(1))]); 
    for n=1:length(IDme), [tmp1, tmp2] = ttest2(res_M1mem(n,:),b2(n,:),0.01); 
        s_vals(n) = tmp1; p_vals(n) = tmp2; end;
    plot(dx(logical(s_vals)),55*ones(size(dx(logical(s_vals)))),'*k')
    [p2, k2] = polyfit(dx(1:end),M2(1:length(IDme)),1); [p1, k1] = polyfit(dx(1:end),M1mem(1:end)',1);
    Z=(p1(1)-p2(1)) / sqrt(k1.normr^2 + k2.normr^2);
    p_slope(6) = normcdf(-abs(Z),0,1);
    hold off; 
    set(gca,'XTick',0:4,'Ytick',45:2:55); title('AB2'); 
    axis([-1 5 45 55]);xlabel('\Deltax');   
    y  = [clm((b2(2:4,:))')'        clm(res_M1mem(2:end,:)')']';
    g1 = [1*ones(1,12)             2*ones(1,12)]; 
    g2 = [clm(((1:3)'*ones(1,4))')' clm(((1:3)'*ones(1,4))')']; 
    group={g1 g2}; if flag_do_stats, [p,table,stats,terms] = anovan(y,group,'model',[1 0;0 1;1 1]); end;
    P(6,:)=p;

subplot(2,4,7)
plot(0,0,'ok','markerfacecolor','k','linewidth',2);hold on;
plot(0,0,'or','markerfacecolor','w','linewidth',2);
plot(0,0,'w');  
plot(0,0,'<',0,0,'>','color',[.5 .5 .5],'markersize',5,'linewidth',2,'markerfacecolor',grey); hold off; xlim([3 4]); axis off;
legend('MLT dual masker', 'MLT most effective','single masker','MLT single basal','MLT single apical','location','North');

set(gcf,'position',[70 10 1050 600])
% print(gcf, '-dpng', '-r200', 'Fig_5_PTC.png');

%% analysis of variance
tmp1=clm(D(2:end,:)); tmp2=clm(S(2:end,:)); 
x=tmp1(tmp1~=0); y=tmp2(tmp2~=0);
g1=[1 1 2 2 2 2 3 3 3 3 4 4 5 5 6 6 6];
g2=[1:2 1:4 1:4 1:2 1:2 1:3];
% aoctool(x,y,[g1' g2']) 
group={g1' g2'}; [p,table,stats,terms] = anovan(y,group);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data from Dingemanse
clear all; clc;
figure(235)


dx=0:5; x=0:0.01:5; xs=0:0.01:3;


%S1
subplot(3,2,1);
S1du=[150.5 133.5 116.5 77.5 43.5 40.5];
S1sa=[118.5 106.5 79.5 51.5];
S1sb=[118.5 94.5 75.5 47.5];
plot(dx,S1du,'-ok',dx(1:4),S1sa,'-<k',dx(1:4),S1sb,'-sk'); hold on;
ydu=interp1(0:5,S1du,x); ysa=interp1(0:3,S1sa,xs); ysb=interp1(0:3,S1sb,xs); yms=interp1(0:3,max(S1sb,S1sa),xs);  
HWx(1,:)=[xs(argmin(abs(yms-S1sa(1)/2))) x(argmin(abs(ydu-S1du(1)/2))) xs(argmin(abs(ysa-S1sa(1)/2))) xs(argmin(abs(ysb-S1sb(1)/2))) x(argmin(abs(ydu-S1du(2)/2))) xs(argmin(abs(ysa-S1sa(2)/2))) xs(argmin(abs(ysb-S1sb(2)/2)))];
HWy(1,:)=[yms(argmin(abs(yms-S1sa(1)/2))) ydu(argmin(abs(ydu-S1du(1)/2))) ysa(argmin(abs(ysa-S1sa(1)/2))) ysb(argmin(abs(ysb-S1sb(1)/2))) ydu(argmin(abs(ydu-S1du(2)/2))) ysa(argmin(abs(ysa-S1sa(2)/2))) ysb(argmin(abs(ysb-S1sb(2)/2)))];
plot([xs(argmin(abs(yms-S1sa(1)/2))) xs(argmin(abs(yms-S1sa(1)/2)))],[yms(argmin(abs(yms-S1sa(1)/2)))-15 yms(argmin(abs(yms-S1sa(1)/2)))+15],'g','linewidth',3);
plot([x(argmin(abs(ydu-S1du(1)/2)))  x(argmin(abs(ydu-S1du(1)/2)))], [ydu(argmin(abs(ydu-S1du(1)/2)))-10 ydu(argmin(abs(ydu-S1du(1)/2)))+10],'k');
plot([xs(argmin(abs(ysa-S1sa(1)/2))) xs(argmin(abs(ysa-S1sa(1)/2)))],[ysa(argmin(abs(ysa-S1sa(1)/2)))-10 ysa(argmin(abs(ysa-S1sa(1)/2)))+10],'k');
plot([xs(argmin(abs(ysb-S1sb(1)/2))) xs(argmin(abs(ysb-S1sb(1)/2)))],[ysb(argmin(abs(ysb-S1sb(1)/2)))-10 ysb(argmin(abs(ysb-S1sb(1)/2)))+10],'k');
plot([x(argmin(abs(ydu-S1du(2)/2)))  x(argmin(abs(ydu-S1du(2)/2)))], [ydu(argmin(abs(ydu-S1du(2)/2)))-10 ydu(argmin(abs(ydu-S1du(2)/2)))+10],'r');
plot([xs(argmin(abs(ysa-S1sa(2)/2))) xs(argmin(abs(ysa-S1sa(2)/2)))],[ysa(argmin(abs(ysa-S1sa(2)/2)))-10 ysa(argmin(abs(ysa-S1sa(2)/2)))+10],'r');
plot([xs(argmin(abs(ysb-S1sb(2)/2))) xs(argmin(abs(ysb-S1sb(2)/2)))],[ysb(argmin(abs(ysb-S1sb(2)/2)))-10 ysb(argmin(abs(ysb-S1sb(2)/2)))+10],'r');
hold off; title('S1');  xlabel('electrodes (and widths)'); ylabel('T shifts [uA]'); 
axis([-.5 5.5 0 180]);
%S2
subplot(3,2,2);
S2du=[140.9 136.9 76.9 58.3 44.9 39.5];
S2sa=[100.9 64.4 50.3 36.9];
S2sb=[100.9 57.8 26.6 22.9];
plot(dx,S2du,'-ok',dx(1:4),S2sa,'-<k',dx(1:4),S2sb,'-sk');hold on;  
ydu=interp1(0:5,S2du,x); ysa=interp1(0:3,S2sa,xs); ysb=interp1(0:3,S2sb,xs); yms=interp1(0:3,max(S2sb,S2sa),xs); 
HWx(2,:)=[xs(argmin(abs(yms-S2sa(1)/2))) x(argmin(abs(ydu-S2du(1)/2))) xs(argmin(abs(ysa-S2sa(1)/2))) xs(argmin(abs(ysb-S2sb(1)/2))) x(argmin(abs(ydu-S2du(2)/2))) xs(argmin(abs(ysa-S2sa(2)/2))) xs(argmin(abs(ysb-S2sb(2)/2)))];
HWy(2,:)=[yms(argmin(abs(yms-S2sa(1)/2))) ydu(argmin(abs(ydu-S2du(1)/2))) ysa(argmin(abs(ysa-S2sa(1)/2))) ysb(argmin(abs(ysb-S2sb(1)/2))) ydu(argmin(abs(ydu-S2du(2)/2))) ysa(argmin(abs(ysa-S2sa(2)/2))) ysb(argmin(abs(ysb-S2sb(2)/2)))];
plot([xs(argmin(abs(yms-S2sa(1)/2))) xs(argmin(abs(yms-S2sa(1)/2)))],[yms(argmin(abs(yms-S2sa(1)/2)))-15 yms(argmin(abs(yms-S2sa(1)/2)))+15],'g','linewidth',3);
plot([x(argmin(abs(ydu-S2du(1)/2))) x(argmin(abs(ydu-S2du(1)/2)))],  [ydu(argmin(abs(ydu-S2du(1)/2)))-10 ydu(argmin(abs(ydu-S2du(1)/2)))+10],'k');
plot([xs(argmin(abs(ysa-S2sa(1)/2))) xs(argmin(abs(ysa-S2sa(1)/2)))],[ysa(argmin(abs(ysa-S2sa(1)/2)))-10 ysa(argmin(abs(ysa-S2sa(1)/2)))+10],'k');
plot([xs(argmin(abs(ysb-S2sb(1)/2))) xs(argmin(abs(ysb-S2sb(1)/2)))],[ysb(argmin(abs(ysb-S2sb(1)/2)))-10 ysb(argmin(abs(ysb-S2sb(1)/2)))+10],'k');
plot([x(argmin(abs(ydu-S2du(2)/2)))  x(argmin(abs(ydu-S2du(2)/2)))], [ydu(argmin(abs(ydu-S2du(2)/2)))-10 ydu(argmin(abs(ydu-S2du(2)/2)))+10],'r');
plot([xs(argmin(abs(ysa-S2sa(2)/2))) xs(argmin(abs(ysa-S2sa(2)/2)))],[ysa(argmin(abs(ysa-S2sa(2)/2)))-10 ysa(argmin(abs(ysa-S2sa(2)/2)))+10],'r');
plot([xs(argmin(abs(ysb-S2sb(2)/2))) xs(argmin(abs(ysb-S2sb(2)/2)))],[ysb(argmin(abs(ysb-S2sb(2)/2)))-10 ysb(argmin(abs(ysb-S2sb(2)/2)))+10],'r');
hold off;title('S2');  xlabel('electrodes (and widths)'); ylabel('T shifts [uA]'); 
axis([-.5 5.5 0 170]);
%S3
subplot(3,2,3);
S3du=[121.5 104.5 73.1 55.5 44.5 39.5];
S3sa=[68.5 46.5 46.5 38.5];
S3sb=[68.5 57.5 41.5 10.5];
plot(dx,S3du,'-ok',dx(1:4),S3sa,'-<k',dx(1:4),S3sb,'-sk'); hold on;
ydu=interp1(0:5,S3du,x); ysa=interp1(0:3,S3sa,xs); ysb=interp1(0:3,S3sb,xs); yms=interp1(0:3,max(S3sb,S3sa),xs);  
HWx(3,:)=[xs(argmin(abs(yms-S3sa(1)/2))) x(argmin(abs(ydu-S3du(1)/2))) xs(argmin(abs(ysa-S3sa(1)/2))) xs(argmin(abs(ysb-S3sb(1)/2))) x(argmin(abs(ydu-S3du(2)/2))) xs(argmin(abs(ysa-S3sa(2)/2))) xs(argmin(abs(ysb-S3sb(2)/2)))];
HWy(3,:)=[yms(argmin(abs(yms-S3sa(1)/2))) ydu(argmin(abs(ydu-S3du(1)/2))) ysa(argmin(abs(ysa-S3sa(1)/2))) ysb(argmin(abs(ysb-S3sb(1)/2))) ydu(argmin(abs(ydu-S3du(2)/2))) ysa(argmin(abs(ysa-S3sa(2)/2))) ysb(argmin(abs(ysb-S3sb(2)/2)))];
plot([xs(argmin(abs(yms-S3sa(1)/2))) xs(argmin(abs(yms-S3sa(1)/2)))],[yms(argmin(abs(yms-S3sa(1)/2)))-15 yms(argmin(abs(yms-S3sa(1)/2)))+15],'g','linewidth',3);
plot([x(argmin(abs(ydu-S3du(1)/2))) x(argmin(abs(ydu-S3du(1)/2)))],  [ydu(argmin(abs(ydu-S3du(1)/2)))-10 ydu(argmin(abs(ydu-S3du(1)/2)))+10],'k');
plot([xs(argmin(abs(ysa-S3sa(1)/2))) xs(argmin(abs(ysa-S3sa(1)/2)))],[ysa(argmin(abs(ysa-S3sa(1)/2)))-10 ysa(argmin(abs(ysa-S3sa(1)/2)))+10],'k');
plot([xs(argmin(abs(ysb-S3sb(1)/2))) xs(argmin(abs(ysb-S3sb(1)/2)))],[ysb(argmin(abs(ysb-S3sb(1)/2)))-10 ysb(argmin(abs(ysb-S3sb(1)/2)))+10],'k');
plot([x(argmin(abs(ydu-S3du(2)/2)))  x(argmin(abs(ydu-S3du(2)/2)))], [ydu(argmin(abs(ydu-S3du(2)/2)))-10 ydu(argmin(abs(ydu-S3du(2)/2)))+10],'r');
plot([xs(argmin(abs(ysa-S3sa(2)/2))) xs(argmin(abs(ysa-S3sa(2)/2)))],[ysa(argmin(abs(ysa-S3sa(2)/2)))-10 ysa(argmin(abs(ysa-S3sa(2)/2)))+10],'r');
plot([xs(argmin(abs(ysb-S3sb(2)/2))) xs(argmin(abs(ysb-S3sb(2)/2)))],[ysb(argmin(abs(ysb-S3sb(2)/2)))-10 ysb(argmin(abs(ysb-S3sb(2)/2)))+10],'r');
hold off; title('S3');  xlabel('electrodes (and widths)'); ylabel('T shifts [uA]'); 
axis([-.5 5.5 0 150]);
%S4
%S5
subplot(3,2,5);
S5du=[197.4 185.5 92.5 31.6 NaN 13.8];
S5sa=[107.5 96.5 42.5 12.0];
S5sb=[107.5 59.5 22.5 0.5];
plot(dx,S5du,'-ok',dx(1:4),S5sa,'-<k',dx(1:4),S5sb,'-sk'); hold on;  
ydu=interp1(0:5,S5du,x); ysa=interp1(0:3,S5sa,xs); ysb=interp1(0:3,S5sb,xs); yms=interp1(0:3,max(S5sb,S5sa),xs); 
HWx(4,:)=[xs(argmin(abs(yms-S5sa(1)/2))) x(argmin(abs(ydu-S5du(1)/2))) xs(argmin(abs(ysa-S5sa(1)/2))) xs(argmin(abs(ysb-S5sb(1)/2))) x(argmin(abs(ydu-S5du(2)/2))) xs(argmin(abs(ysa-S5sa(2)/2))) xs(argmin(abs(ysb-S5sb(2)/2)))];
HWy(4,:)=[yms(argmin(abs(yms-S5sa(1)/2))) ydu(argmin(abs(ydu-S5du(1)/2))) ysa(argmin(abs(ysa-S5sa(1)/2))) ysb(argmin(abs(ysb-S5sb(1)/2))) ydu(argmin(abs(ydu-S5du(2)/2))) ysa(argmin(abs(ysa-S5sa(2)/2))) ysb(argmin(abs(ysb-S5sb(2)/2)))];
plot([xs(argmin(abs(yms-S5sa(1)/2))) xs(argmin(abs(yms-S5sa(1)/2)))],[yms(argmin(abs(yms-S5sa(1)/2)))-15 yms(argmin(abs(yms-S5sa(1)/2)))+15],'g','linewidth',3);
plot([x(argmin(abs(ydu-S5du(1)/2))) x(argmin(abs(ydu-S5du(1)/2)))],  [ydu(argmin(abs(ydu-S5du(1)/2)))-10 ydu(argmin(abs(ydu-S5du(1)/2)))+10],'k');
plot([xs(argmin(abs(ysa-S5sa(1)/2))) xs(argmin(abs(ysa-S5sa(1)/2)))],[ysa(argmin(abs(ysa-S5sa(1)/2)))-10 ysa(argmin(abs(ysa-S5sa(1)/2)))+10],'k');
plot([xs(argmin(abs(ysb-S5sb(1)/2))) xs(argmin(abs(ysb-S5sb(1)/2)))],[ysb(argmin(abs(ysb-S5sb(1)/2)))-10 ysb(argmin(abs(ysb-S5sb(1)/2)))+10],'k');
plot([x(argmin(abs(ydu-S5du(2)/2)))  x(argmin(abs(ydu-S5du(2)/2)))], [ydu(argmin(abs(ydu-S5du(2)/2)))-10 ydu(argmin(abs(ydu-S5du(2)/2)))+10],'r');
plot([xs(argmin(abs(ysa-S5sa(2)/2))) xs(argmin(abs(ysa-S5sa(2)/2)))],[ysa(argmin(abs(ysa-S5sa(2)/2)))-10 ysa(argmin(abs(ysa-S5sa(2)/2)))+10],'r');
plot([xs(argmin(abs(ysb-S5sb(2)/2))) xs(argmin(abs(ysb-S5sb(2)/2)))],[ysb(argmin(abs(ysb-S5sb(2)/2)))-10 ysb(argmin(abs(ysb-S5sb(2)/2)))+10],'r');
hold off; title('S5');  xlabel('electrodes (and widths)'); ylabel('T shifts [uA]'); 
axis([-.5 5.5 0 250]);
%S6
subplot(3,2,6); title('S6');
S6du=[196.5 172.5 120.5 89.5 62.5 62.5];
S6sa=[142.5 91.5 51.5 34.5];
S6sb=[142.5 117.5 79.5 62.5];
plot(dx,S6du,'-ok',dx(1:4),S6sa,'-<k',dx(1:4),S6sb,'-sk'); hold on; 
ydu=interp1(0:5,S6du,x); ysa=interp1(0:3,S6sa,xs); ysb=interp1(0:3,S6sb,xs); yms=interp1(0:3,max(S6sb,S6sa),xs);
HWx(5,:)=[xs(argmin(abs(yms-S6sa(1)/2))) x(argmin(abs(ydu-S6du(1)/2))) xs(argmin(abs(ysa-S6sa(1)/2))) xs(argmin(abs(ysb-S6sb(1)/2))) x(argmin(abs(ydu-S6du(2)/2))) xs(argmin(abs(ysa-S6sa(2)/2))) xs(argmin(abs(ysb-S6sb(2)/2)))];
HWy(5,:)=[yms(argmin(abs(yms-S6sa(1)/2))) ydu(argmin(abs(ydu-S6du(1)/2))) ysa(argmin(abs(ysa-S6sa(1)/2))) ysb(argmin(abs(ysb-S6sb(1)/2))) ydu(argmin(abs(ydu-S6du(2)/2))) ysa(argmin(abs(ysa-S6sa(2)/2))) ysb(argmin(abs(ysb-S6sb(2)/2)))];
plot([xs(argmin(abs(yms-S6sa(1)/2))) xs(argmin(abs(yms-S6sa(1)/2)))],[yms(argmin(abs(yms-S6sa(1)/2)))-15 yms(argmin(abs(yms-S6sa(1)/2)))+15],'g','linewidth',3);
plot([ x(argmin(abs(ydu-S6du(1)/2)))  x(argmin(abs(ydu-S6du(1)/2)))],[ydu(argmin(abs(ydu-S6du(1)/2)))-10 ydu(argmin(abs(ydu-S6du(1)/2)))+10],'k','linewidth',1);
plot([xs(argmin(abs(ysa-S6sa(1)/2))) xs(argmin(abs(ysa-S6sa(1)/2)))],[ysa(argmin(abs(ysa-S6sa(1)/2)))-10 ysa(argmin(abs(ysa-S6sa(1)/2)))+10],'k');
plot([xs(argmin(abs(ysb-S6sb(1)/2))) xs(argmin(abs(ysb-S6sb(1)/2)))],[ysb(argmin(abs(ysb-S6sb(1)/2)))-10 ysb(argmin(abs(ysb-S6sb(1)/2)))+10],'k');
plot([ x(argmin(abs(ydu-S6du(2)/2)))  x(argmin(abs(ydu-S6du(2)/2)))],[ydu(argmin(abs(ydu-S6du(2)/2)))-10 ydu(argmin(abs(ydu-S6du(2)/2)))+10],'r');
plot([xs(argmin(abs(ysa-S6sa(2)/2))) xs(argmin(abs(ysa-S6sa(2)/2)))],[ysa(argmin(abs(ysa-S6sa(2)/2)))-10 ysa(argmin(abs(ysa-S6sa(2)/2)))+10],'r');
plot([xs(argmin(abs(ysb-S6sb(2)/2))) xs(argmin(abs(ysb-S6sb(2)/2)))],[ysb(argmin(abs(ysb-S6sb(2)/2)))-10 ysb(argmin(abs(ysb-S6sb(2)/2)))+10],'r');
hold off; xlabel('electrodes (and widths)'); ylabel('T shifts [uA]'); 
title('S6'); axis([-.5 5.5 0 250]); 

subplot(3,2,4);
plot(-1,-1,'ok',-1,-1,'<k',-1,-1,'sk',-1,-1,'-k',-1,-1,'-r',-1,-1,'-g',-1,-1,'w');
legend('dual', 'apical', 'basal','HW as in Dingemanse', 'HW excluding point Dx=0','HW for most effective single', 'masker (Dx=0 included)','location','North');
axis off;

set(gcf,'position',[70 10 600 1050])


figure(236)
nn=[1:3 5 6];
subplot(211), for n=1:5, plot(1,HWx(n,2),'ok',...
                              2,HWx(n,3),'<k',...
                              3,HWx(n,4),'sk',...
                              1:3,HWx(n,2:4),'--k'); 
                              text(.8,HWx(n,2),['S',num2str(nn(n))]);
                              hold on; end; hold off;
                set(gca,'xtick',1:3,'xticklabel',{'dual' 'apical' 'basal'}); xlim([.5 3.5]);ylabel('half widths'); 
subplot(212), for n=1:5, plot(n,HWx(n,2),'ok',n,HWx(n,3),'<k',n,HWx(n,4),'sk'); hold on; end; hold off; ylabel('half widths'); 
    xlim([.5 5.5]); set(gca,'xtick',1:5,'xticklabel',{'S1' 'S2' 'S3' 'S5' 'S6'}); legend('dual masker','single apical','single basal');
set(gcf,'position',[680 10 550 1050])
%%% estimate HWs and statistical significance
clear res;
[h, p] = ttest(HWx(:,3),HWx(:,4)); res(1) = p; % apical Vs basal
[h, p] = ttest(HWx(:,3),HWx(:,2)); res(2) = p; % apical Vs dual
[h, p] = ttest(HWx(:,4),HWx(:,2)); res(3) = p; % basal Vs dual
[h, p] = ttest(nanmin(HWx(:,[3 4])')',HWx(:,2)); res(4) = p; % least effective Vs dual
[h, p] = ttest(nanmax(HWx(:,[3 4])')',HWx(:,2)); res(5) = p; % most effective Vs dual
[h, p] = ttest(sum(HWx(:,[3 4])')',2*HWx(:,2)); res(7) = p; % most effective Vs dual
[h, p] = ttest(HWx(:,3+3),HWx(:,4+3)); res2(1) = p; % apical Vs basal
[h, p] = ttest(HWx(:,3+3),HWx(:,2+3)); res2(2) = p; % apical Vs dual
[h, p] = ttest(HWx(:,4+3),HWx(:,2+3)); res2(3) = p; % basal Vs dual
[h, p] = ttest(nanmin(HWx(:,[3+3 4+3])')',HWx(:,2+3)); res2(4) = p; % least effective Vs dual
[h, p] = ttest(nanmax(HWx(:,[3+3 4+3])')',HWx(:,2+3)); res2(5) = p; % most effective Vs dual
[h, p] = ttest(sum(HWx(:,[3+3 4+3])')',2*HWx(:,2+3)); res2(7) = p; % most effective Vs dual
HWx2=HWx; HWx2(3,3)=HWx(3,4); % use S3 basal like in paper
[h, p] = ttest(nanmax(HWx2(:,[3 4])')',HWx2(:,2)); res(6) = p; % most effective Vs dual
[h, p] = ttest(sum(HWx2(:,[3 4])')',2*HWx2(:,2)); res(8) = p; % most effective Vs dual
[h, p] = ttest(nanmax(HWx2(:,[3+3 4+3])')',HWx2(:,2+3)); res2(6) = p; % most effective Vs dual
[h, p] = ttest(sum(HWx2(:,[3+3 4+3])')',2*HWx2(:,2+3)); res2(8) = p; % most effective Vs dual
round(1000*res)/1000
round(1000*res2)/1000

%% Dingemanse in dB

clear all; clc;
figure(235)


dx=0:5; x=0:0.01:5; xs=0:0.01:3;
Plev=[245 297 371 276 218]; % probe levels, from Table 2
S1du=20*log10([150.5 133.5 116.5 77.5 43.5 40.5]+Plev(1))-20*log10(Plev(1));
S1sa=20*log10([118.5 106.5 79.5 51.5]+Plev(1))-20*log10(Plev(1));
S1sb=20*log10([118.5 94.5 75.5 47.5]+Plev(1))-20*log10(Plev(1));
S2du=20*log10([140.9 136.9 76.9 58.3 44.9 39.5]+Plev(2))-20*log10(Plev(2));
S2sa=20*log10([100.9 64.4 50.3 36.9]+Plev(2))-20*log10(Plev(2));
S2sb=20*log10([100.9 57.8 26.6 22.9]+Plev(2))-20*log10(Plev(2));
S3du=20*log10([121.5 104.5 73.1 55.5 44.5 39.5]+Plev(3))-20*log10(Plev(3));
S3sa=20*log10([68.5 46.5 46.5 38.5]+Plev(3))-20*log10(Plev(3));
S3sb=20*log10([68.5 57.5 41.5 10.5]+Plev(3))-20*log10(Plev(3));
S5du=20*log10([197.4 185.5 92.5 31.6 NaN 13.8]+Plev(4))-20*log10(Plev(4));
S5sa=20*log10([107.5 96.5 42.5 12.0]+Plev(4))-20*log10(Plev(4));
S5sb=20*log10([107.5 59.5 22.5 0.5]+Plev(4))-20*log10(Plev(4));
S6du=20*log10([196.5 172.5 120.5 89.5 62.5 62.5]+Plev(5))-20*log10(Plev(5));
S6sa=20*log10([142.5 91.5 51.5 34.5]+Plev(5))-20*log10(Plev(5));
S6sb=20*log10([142.5 117.5 79.5 62.5]+Plev(5))-20*log10(Plev(5));


%S1
subplot(3,2,1);
plot(dx,S1du,'-ok',dx(1:4),S1sa,'-<k',dx(1:4),S1sb,'-sk'); hold on;
ydu=interp1(0:5,S1du,x); ysa=interp1(0:3,S1sa,xs); ysb=interp1(0:3,S1sb,xs); yms=interp1(0:3,max(S1sb,S1sa),xs);  
HWx(1,:)=[xs(argmin(abs(yms-S1sa(1)/2))) x(argmin(abs(ydu-S1du(1)/2))) xs(argmin(abs(ysa-S1sa(1)/2))) xs(argmin(abs(ysb-S1sb(1)/2))) x(argmin(abs(ydu-S1du(2)/2))) xs(argmin(abs(ysa-S1sa(2)/2))) xs(argmin(abs(ysb-S1sb(2)/2)))];
HWy(1,:)=[yms(argmin(abs(yms-S1sa(1)/2))) ydu(argmin(abs(ydu-S1du(1)/2))) ysa(argmin(abs(ysa-S1sa(1)/2))) ysb(argmin(abs(ysb-S1sb(1)/2))) ydu(argmin(abs(ydu-S1du(2)/2))) ysa(argmin(abs(ysa-S1sa(2)/2))) ysb(argmin(abs(ysb-S1sb(2)/2)))];
plot([xs(argmin(abs(yms-S1sa(1)/2))) xs(argmin(abs(yms-S1sa(1)/2)))],[yms(argmin(abs(yms-S1sa(1)/2)))-.15 yms(argmin(abs(yms-S1sa(1)/2)))+.15],'g','linewidth',3);
plot([x(argmin(abs(ydu-S1du(1)/2)))  x(argmin(abs(ydu-S1du(1)/2)))], [ydu(argmin(abs(ydu-S1du(1)/2)))-.10 ydu(argmin(abs(ydu-S1du(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysa-S1sa(1)/2))) xs(argmin(abs(ysa-S1sa(1)/2)))],[ysa(argmin(abs(ysa-S1sa(1)/2)))-.10 ysa(argmin(abs(ysa-S1sa(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysb-S1sb(1)/2))) xs(argmin(abs(ysb-S1sb(1)/2)))],[ysb(argmin(abs(ysb-S1sb(1)/2)))-.10 ysb(argmin(abs(ysb-S1sb(1)/2)))+.10],'k');
plot([x(argmin(abs(ydu-S1du(2)/2)))  x(argmin(abs(ydu-S1du(2)/2)))], [ydu(argmin(abs(ydu-S1du(2)/2)))-.10 ydu(argmin(abs(ydu-S1du(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysa-S1sa(2)/2))) xs(argmin(abs(ysa-S1sa(2)/2)))],[ysa(argmin(abs(ysa-S1sa(2)/2)))-.10 ysa(argmin(abs(ysa-S1sa(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysb-S1sb(2)/2))) xs(argmin(abs(ysb-S1sb(2)/2)))],[ysb(argmin(abs(ysb-S1sb(2)/2)))-.10 ysb(argmin(abs(ysb-S1sb(2)/2)))+.10],'r');
hold off; title('S1');  xlabel('electrodes (and widths)'); ylabel('T shifts [db]'); 
% axis([-.5 5.5 0 180]);
%S2
subplot(3,2,2);
plot(dx,S2du,'-ok',dx(1:4),S2sa,'-<k',dx(1:4),S2sb,'-sk');hold on;  
ydu=interp1(0:5,S2du,x); ysa=interp1(0:3,S2sa,xs); ysb=interp1(0:3,S2sb,xs); yms=interp1(0:3,max(S2sb,S2sa),xs); 
HWx(2,:)=[xs(argmin(abs(yms-S2sa(1)/2))) x(argmin(abs(ydu-S2du(1)/2))) xs(argmin(abs(ysa-S2sa(1)/2))) xs(argmin(abs(ysb-S2sb(1)/2))) x(argmin(abs(ydu-S2du(2)/2))) xs(argmin(abs(ysa-S2sa(2)/2))) xs(argmin(abs(ysb-S2sb(2)/2)))];
HWy(2,:)=[yms(argmin(abs(yms-S2sa(1)/2))) ydu(argmin(abs(ydu-S2du(1)/2))) ysa(argmin(abs(ysa-S2sa(1)/2))) ysb(argmin(abs(ysb-S2sb(1)/2))) ydu(argmin(abs(ydu-S2du(2)/2))) ysa(argmin(abs(ysa-S2sa(2)/2))) ysb(argmin(abs(ysb-S2sb(2)/2)))];
plot([xs(argmin(abs(yms-S2sa(1)/2))) xs(argmin(abs(yms-S2sa(1)/2)))],[yms(argmin(abs(yms-S2sa(1)/2)))-.15 yms(argmin(abs(yms-S2sa(1)/2)))+.15],'g','linewidth',3);
plot([x(argmin(abs(ydu-S2du(1)/2))) x(argmin(abs(ydu-S2du(1)/2)))],  [ydu(argmin(abs(ydu-S2du(1)/2)))-.10 ydu(argmin(abs(ydu-S2du(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysa-S2sa(1)/2))) xs(argmin(abs(ysa-S2sa(1)/2)))],[ysa(argmin(abs(ysa-S2sa(1)/2)))-.10 ysa(argmin(abs(ysa-S2sa(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysb-S2sb(1)/2))) xs(argmin(abs(ysb-S2sb(1)/2)))],[ysb(argmin(abs(ysb-S2sb(1)/2)))-.10 ysb(argmin(abs(ysb-S2sb(1)/2)))+.10],'k');
plot([x(argmin(abs(ydu-S2du(2)/2)))  x(argmin(abs(ydu-S2du(2)/2)))], [ydu(argmin(abs(ydu-S2du(2)/2)))-.10 ydu(argmin(abs(ydu-S2du(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysa-S2sa(2)/2))) xs(argmin(abs(ysa-S2sa(2)/2)))],[ysa(argmin(abs(ysa-S2sa(2)/2)))-.10 ysa(argmin(abs(ysa-S2sa(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysb-S2sb(2)/2))) xs(argmin(abs(ysb-S2sb(2)/2)))],[ysb(argmin(abs(ysb-S2sb(2)/2)))-.10 ysb(argmin(abs(ysb-S2sb(2)/2)))+.10],'r');
hold off;title('S2');  xlabel('electrodes (and widths)'); ylabel('T shifts [db]'); 
% axis([-.5 5.5 0 170]);
%S3
subplot(3,2,3);
plot(dx,S3du,'-ok',dx(1:4),S3sa,'-<k',dx(1:4),S3sb,'-sk'); hold on;
ydu=interp1(0:5,S3du,x); ysa=interp1(0:3,S3sa,xs); ysb=interp1(0:3,S3sb,xs); yms=interp1(0:3,max(S3sb,S3sa),xs);  
HWx(3,:)=[xs(argmin(abs(yms-S3sa(1)/2))) x(argmin(abs(ydu-S3du(1)/2))) xs(argmin(abs(ysa-S3sa(1)/2))) xs(argmin(abs(ysb-S3sb(1)/2))) x(argmin(abs(ydu-S3du(2)/2))) xs(argmin(abs(ysa-S3sa(2)/2))) xs(argmin(abs(ysb-S3sb(2)/2)))];
HWy(3,:)=[yms(argmin(abs(yms-S3sa(1)/2))) ydu(argmin(abs(ydu-S3du(1)/2))) ysa(argmin(abs(ysa-S3sa(1)/2))) ysb(argmin(abs(ysb-S3sb(1)/2))) ydu(argmin(abs(ydu-S3du(2)/2))) ysa(argmin(abs(ysa-S3sa(2)/2))) ysb(argmin(abs(ysb-S3sb(2)/2)))];
plot([xs(argmin(abs(yms-S3sa(1)/2))) xs(argmin(abs(yms-S3sa(1)/2)))],[yms(argmin(abs(yms-S3sa(1)/2)))-.15 yms(argmin(abs(yms-S3sa(1)/2)))+.15],'g','linewidth',3);
plot([x(argmin(abs(ydu-S3du(1)/2))) x(argmin(abs(ydu-S3du(1)/2)))],  [ydu(argmin(abs(ydu-S3du(1)/2)))-.10 ydu(argmin(abs(ydu-S3du(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysa-S3sa(1)/2))) xs(argmin(abs(ysa-S3sa(1)/2)))],[ysa(argmin(abs(ysa-S3sa(1)/2)))-.10 ysa(argmin(abs(ysa-S3sa(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysb-S3sb(1)/2))) xs(argmin(abs(ysb-S3sb(1)/2)))],[ysb(argmin(abs(ysb-S3sb(1)/2)))-.10 ysb(argmin(abs(ysb-S3sb(1)/2)))+.10],'k');
plot([x(argmin(abs(ydu-S3du(2)/2)))  x(argmin(abs(ydu-S3du(2)/2)))], [ydu(argmin(abs(ydu-S3du(2)/2)))-.10 ydu(argmin(abs(ydu-S3du(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysa-S3sa(2)/2))) xs(argmin(abs(ysa-S3sa(2)/2)))],[ysa(argmin(abs(ysa-S3sa(2)/2)))-.10 ysa(argmin(abs(ysa-S3sa(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysb-S3sb(2)/2))) xs(argmin(abs(ysb-S3sb(2)/2)))],[ysb(argmin(abs(ysb-S3sb(2)/2)))-.10 ysb(argmin(abs(ysb-S3sb(2)/2)))+.10],'r');
hold off; title('S3');  xlabel('electrodes (and widths)'); ylabel('T shifts [db]'); 
% axis([-.5 5.5 0 150]);
%S4
%S5
subplot(3,2,5);
plot(dx,S5du,'-ok',dx(1:4),S5sa,'-<k',dx(1:4),S5sb,'-sk'); hold on;  
ydu=interp1(0:5,S5du,x); ysa=interp1(0:3,S5sa,xs); ysb=interp1(0:3,S5sb,xs); yms=interp1(0:3,max(S5sb,S5sa),xs); 
HWx(4,:)=[xs(argmin(abs(yms-S5sa(1)/2))) x(argmin(abs(ydu-S5du(1)/2))) xs(argmin(abs(ysa-S5sa(1)/2))) xs(argmin(abs(ysb-S5sb(1)/2))) x(argmin(abs(ydu-S5du(2)/2))) xs(argmin(abs(ysa-S5sa(2)/2))) xs(argmin(abs(ysb-S5sb(2)/2)))];
HWy(4,:)=[yms(argmin(abs(yms-S5sa(1)/2))) ydu(argmin(abs(ydu-S5du(1)/2))) ysa(argmin(abs(ysa-S5sa(1)/2))) ysb(argmin(abs(ysb-S5sb(1)/2))) ydu(argmin(abs(ydu-S5du(2)/2))) ysa(argmin(abs(ysa-S5sa(2)/2))) ysb(argmin(abs(ysb-S5sb(2)/2)))];
plot([xs(argmin(abs(yms-S5sa(1)/2))) xs(argmin(abs(yms-S5sa(1)/2)))],[yms(argmin(abs(yms-S5sa(1)/2)))-.15 yms(argmin(abs(yms-S5sa(1)/2)))+.15],'g','linewidth',3);
plot([x(argmin(abs(ydu-S5du(1)/2)))   x(argmin(abs(ydu-S5du(1)/2)))],[ydu(argmin(abs(ydu-S5du(1)/2)))-.10 ydu(argmin(abs(ydu-S5du(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysa-S5sa(1)/2))) xs(argmin(abs(ysa-S5sa(1)/2)))],[ysa(argmin(abs(ysa-S5sa(1)/2)))-.10 ysa(argmin(abs(ysa-S5sa(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysb-S5sb(1)/2))) xs(argmin(abs(ysb-S5sb(1)/2)))],[ysb(argmin(abs(ysb-S5sb(1)/2)))-.10 ysb(argmin(abs(ysb-S5sb(1)/2)))+.10],'k');
plot([x(argmin(abs(ydu-S5du(2)/2)))  x(argmin(abs(ydu-S5du(2)/2)))], [ydu(argmin(abs(ydu-S5du(2)/2)))-.10 ydu(argmin(abs(ydu-S5du(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysa-S5sa(2)/2))) xs(argmin(abs(ysa-S5sa(2)/2)))],[ysa(argmin(abs(ysa-S5sa(2)/2)))-.10 ysa(argmin(abs(ysa-S5sa(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysb-S5sb(2)/2))) xs(argmin(abs(ysb-S5sb(2)/2)))],[ysb(argmin(abs(ysb-S5sb(2)/2)))-.10 ysb(argmin(abs(ysb-S5sb(2)/2)))+.10],'r');
hold off; title('S5');  xlabel('electrodes (and widths)'); ylabel('T shifts [db]'); 
% axis([-.5 5.5 0 250]);
%S6
subplot(3,2,6); title('S6');
plot(dx,S6du,'-ok',dx(1:4),S6sa,'-<k',dx(1:4),S6sb,'-sk'); hold on; 
ydu=interp1(0:5,S6du,x); ysa=interp1(0:3,S6sa,xs); ysb=interp1(0:3,S6sb,xs); yms=interp1(0:3,max(S6sb,S6sa),xs);
HWx(5,:)=[xs(argmin(abs(yms-S6sa(1)/2))) x(argmin(abs(ydu-S6du(1)/2))) xs(argmin(abs(ysa-S6sa(1)/2))) xs(argmin(abs(ysb-S6sb(1)/2))) x(argmin(abs(ydu-S6du(2)/2))) xs(argmin(abs(ysa-S6sa(2)/2))) xs(argmin(abs(ysb-S6sb(2)/2)))];
HWy(5,:)=[yms(argmin(abs(yms-S6sa(1)/2))) ydu(argmin(abs(ydu-S6du(1)/2))) ysa(argmin(abs(ysa-S6sa(1)/2))) ysb(argmin(abs(ysb-S6sb(1)/2))) ydu(argmin(abs(ydu-S6du(2)/2))) ysa(argmin(abs(ysa-S6sa(2)/2))) ysb(argmin(abs(ysb-S6sb(2)/2)))];
plot([xs(argmin(abs(yms-S6sa(1)/2))) xs(argmin(abs(yms-S6sa(1)/2)))],[yms(argmin(abs(yms-S6sa(1)/2)))-.15 yms(argmin(abs(yms-S6sa(1)/2)))+.15],'g','linewidth',3);
plot([ x(argmin(abs(ydu-S6du(1)/2)))  x(argmin(abs(ydu-S6du(1)/2)))],[ydu(argmin(abs(ydu-S6du(1)/2)))-.10 ydu(argmin(abs(ydu-S6du(1)/2)))+.10],'k','linewidth',1);
plot([xs(argmin(abs(ysa-S6sa(1)/2))) xs(argmin(abs(ysa-S6sa(1)/2)))],[ysa(argmin(abs(ysa-S6sa(1)/2)))-.10 ysa(argmin(abs(ysa-S6sa(1)/2)))+.10],'k');
plot([xs(argmin(abs(ysb-S6sb(1)/2))) xs(argmin(abs(ysb-S6sb(1)/2)))],[ysb(argmin(abs(ysb-S6sb(1)/2)))-.10 ysb(argmin(abs(ysb-S6sb(1)/2)))+.10],'k');
plot([ x(argmin(abs(ydu-S6du(2)/2)))  x(argmin(abs(ydu-S6du(2)/2)))],[ydu(argmin(abs(ydu-S6du(2)/2)))-.10 ydu(argmin(abs(ydu-S6du(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysa-S6sa(2)/2))) xs(argmin(abs(ysa-S6sa(2)/2)))],[ysa(argmin(abs(ysa-S6sa(2)/2)))-.10 ysa(argmin(abs(ysa-S6sa(2)/2)))+.10],'r');
plot([xs(argmin(abs(ysb-S6sb(2)/2))) xs(argmin(abs(ysb-S6sb(2)/2)))],[ysb(argmin(abs(ysb-S6sb(2)/2)))-.10 ysb(argmin(abs(ysb-S6sb(2)/2)))+.10],'r');
hold off; xlabel('electrodes (and widths)'); ylabel('T shifts [db]'); 
title('S6'); 
% axis([-.5 5.5 0 250]); 

subplot(3,2,4);
plot(-1,-1,'ok',-1,-1,'<k',-1,-1,'sk',-1,-1,'-k',-1,-1,'-r',-1,-1,'-g',-1,-1,'w');
legend('dual', 'apical', 'basal','HW as in Dingemanse', 'HW excluding point Dx=0','HW for most effective single', 'masker (Dx=0 included)','location','North');
axis off;

set(gcf,'position',[70 10 600 1050])


figure(236)
nn=[1:3 5 6];
subplot(211), for n=1:5, plot(1,HWx(n,2),'ok',...
                              2,HWx(n,3),'<k',...
                              3,HWx(n,4),'sk',...
                              1:3,HWx(n,2:4),'--k'); 
                              text(.8,HWx(n,2),['S',num2str(nn(n))]);
                              hold on; end; hold off;
                set(gca,'xtick',1:3,'xticklabel',{'dual' 'apical' 'basal'}); xlim([.5 3.5]);ylabel('half widths'); 
subplot(212), for n=1:5, plot(n,HWx(n,2),'ok',n,HWx(n,3),'<k',n,HWx(n,4),'sk'); hold on; end; hold off; ylabel('half widths'); 
    xlim([.5 5.5]); set(gca,'xtick',1:5,'xticklabel',{'S1' 'S2' 'S3' 'S5' 'S6'}); legend('dual masker','single apical','single basal');
set(gcf,'position',[680 10 550 1050])

%%% estimate HWs and statistical significance
clear res;
[h, p] = ttest(HWx(:,3),HWx(:,4)); res(1) = p; % apical Vs basal
[h, p] = ttest(HWx(:,3),HWx(:,2)); res(2) = p; % apical Vs dual
[h, p] = ttest(HWx(:,4),HWx(:,2)); res(3) = p; % basal Vs dual
[h, p] = ttest(nanmin(HWx(:,[3 4])')',HWx(:,2)); res(4) = p; % least effective Vs dual
[h, p] = ttest(nanmax(HWx(:,[3 4])')',HWx(:,2)); res(5) = p; % most effective Vs dual
[h, p] = ttest(sum(HWx(:,[3 4])')',2*HWx(:,2)); res(7) = p; % most effective Vs dual
[h, p] = ttest(HWx(:,3+3),HWx(:,4+3)); res2(1) = p; % apical Vs basal
[h, p] = ttest(HWx(:,3+3),HWx(:,2+3)); res2(2) = p; % apical Vs dual
[h, p] = ttest(HWx(:,4+3),HWx(:,2+3)); res2(3) = p; % basal Vs dual
[h, p] = ttest(nanmin(HWx(:,[3+3 4+3])')',HWx(:,2+3)); res2(4) = p; % least effective Vs dual
[h, p] = ttest(nanmax(HWx(:,[3+3 4+3])')',HWx(:,2+3)); res2(5) = p; % most effective Vs dual
[h, p] = ttest(sum(HWx(:,[3+3 4+3])')',2*HWx(:,2+3)); res2(7) = p; % most effective Vs dual
HWx(3,3)=HWx(3,4); % use S3 basal like in paper
[h, p] = ttest(nanmax(HWx(:,[3 4])')',HWx(:,2)); res(6) = p; % most effective Vs dual
[h, p] = ttest(sum(HWx(:,[3 4])')',2*HWx(:,2)); res(8) = p; % most effective Vs dual
[h, p] = ttest(nanmax(HWx(:,[3+3 4+3])')',HWx(:,2+3)); res2(6) = p; % most effective Vs dual
[h, p] = ttest(sum(HWx(:,[3+3 4+3])')',2*HWx(:,2+3)); res2(8) = p; % most effective Vs dual
round(1000*res)/1000
round(1000*res2)/1000

