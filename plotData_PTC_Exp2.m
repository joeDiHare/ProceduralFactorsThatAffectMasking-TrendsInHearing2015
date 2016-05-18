%%% Exp. 2 in PTC article
%%% Comparing 2 procedures in AB1
% close all; 
clc; clear all;

Tp([6 8]) = 20*log10([320.0 330.0]);
Cp([6 8]) = 20*log10([273.7 560.0]); 

Lp([6 8]) = 20*log10([345.0 368.0]);

T2(8:12) = 20*log10([240.0 233.7 245.1 261.8 269.4]);
C2(8:12) = 20*log10([390.0 390.0 390.0 405.0 440.0]); %C8,8?

T1(5:11) = 20*log10([269.8 252.9 248.6 294.6 311.5 323.5 351.9]);
C1(5:11) = 20*log10([450.0 450.0 440.0 480.0 490.0 530.0 510.0]); 

M2_8(8:11) = 20*log10([314.7 290.0 313.8 376.0]);
E2_8(8:11) = [0.044 0.192 0.137 0.000];

M2_6(8:11) = 20*log10([290.2 307.4 356.2 379.6]);
E2_6(8:11) = [0.295 0.258 0.407 0.398];

X2_8(8:11) = 20*log10([408.3 400.7 385.4 351.9]);
S2_8(8:11) = [0.254 0.341 0.185 0.133];

M1_8(6:9) = 20*log10([438.9 418.9 469.0 512.7]);
E1_8(6:9) = [0.043 0.540 0.258 0.088];


close(figure(91)); figure(91);

%Subj1 AnthonyM single on E8
subplot(121)
plot(0:3,T1(8:11),'>b',0:3,C1(8:11),'^b',...
     0:3,T1(8:-1:5),'<b',0:3,C1(8:-1:5),'vb','markersize',5,'color',grey,'markerfacecolor',grey); hold on;
% line([0 4],[Tp(8) Tp(8)],'linestyle','--','color',[.5 .5 .5]);
% line([0 4],[Cp(8) Cp(8)],'linestyle','--','color',[.5 .5 .5]); 
errorbar(0:3,X2_8(8:11),S2_8(8:11),'o-k','linewidth',2.5,'markerfacecolor','k'); hold off;
% text(-.5,54.5,'E8','edgecolor','k');
% text(3,Cp(8)-.2,'Cprobe','color',[.5 .5 .5]);
% text(3,Tp(8)+.15,'Tprobe','color',[.5 .5 .5]);
set(gca,'XTick',0:3,'xlim',[-1 4],'ylim',[47.0 55.5]); 
title('Method B'); ylabel('Probe level at threshold [dB]'); xlabel('\Deltax'); 
% h_legend=legend('T1 basal','C1 basal','T1 apical','C1 apical','location','southeast'); 
% set(h_legend,'FontSize',8,'Position',[0.8679 0.1317 0.1275 0.1783]);

%Subj1 AnthonyM probe on E8
subplot(122)
% errorbar(0:3,M2_6(8:11),E2_6(8:11),'o-m','linewidth',2,'color',[.5 .5 .5]);
errorbar(0:3,M2_8(8:11),E2_8(8:11),'o-k','markerfacecolor','k','linewidth',2.5); hold on;
plot(0:3,T2(8:11),'vb',0:3,C2(8:11),'^b','markersize',5,'color',grey,'markerfacecolor',grey); 
% line([0 4],[Lp(6) Lp(6)],'linestyle','--','color',[.5 .5 .5]);
% line([0 4],[Lp(8) Lp(8)],'linestyle','--','color','k'); hold off;
% text(4-.1,52.8,'*','color','k','fontsize',18);
% text(-.5,54.5,'E8','edgecolor','k'); %text(-.5,53.8,'E6','edgecolor',[.5 .5 .5]);
set(gca,'XTick',0:3,'xlim',[-1 4],'ylim',[47.0 55.5]); 
title('Method C'); ylabel('Maskers level at threshold [dB]'); xlabel('\Deltax'); 

set(gcf,'position',[100 10 800 350]);
print(gcf, '-dpng', '-r200', 'Fig_6_PTC.png');


%%
figure(102)
errorbar((0:3)+0.01,-X2_8(8:11)+X2_8(8),S2_8(8:11),'o-','linewidth',2,'color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5]);hold on;
errorbar((0:3)-0.01,M2_8(8:11)-M2_8(8),E2_8(8:11),'s-k','linewidth',2,'markerfacecolor','k'); hold off;
legend('method B','method C','location','northwest');
set(gca,'XTick',0:4,'xlim',[-1 5],'ylim',[-1 2]); 
ylabel('Masking Re:\Deltax=0 [dB]');  
xlabel({'\Deltax'}); 
set(gcf,'position',[100 10 400 350]);


%% AB SUbj1 (AM) 4 plots
% close all; 
clc; clear all;

Tp([6 8]) = 20*log10([320.0 330.0]);
Cp([6 8]) = 20*log10([273.7 560.0]); 

Lp([6 8]) = 20*log10([345.0 368.0]);

T2(8:12) = 20*log10([240.0 233.7 245.1 261.8 269.4]);
C2(8:12) = 20*log10([390.0 390.0 390.0 405.0 440.0]); %C8,8?

T1(5:11) = 20*log10([269.8 252.9 248.6 294.6 311.5 323.5 351.9]);
C1(5:11) = 20*log10([450.0 450.0 440.0 480.0 490.0 530.0 510.0]); 

M2_8(8:11) = 20*log10([314.7 290.0 313.8 376.0]);
E2_8(8:11) = [0.044 0.192 0.137 0.000];

M2_6(8:11) = 20*log10([290.2 307.4 356.2 379.6]);
E2_6(8:11) = [0.295 0.258 0.407 0.398];

X2_8(8:11) = 20*log10([408.3 400.7 385.4 351.9]);
S2_8(8:11) = [0.254 0.341 0.185 0.133];

M1_8(6:9) = 20*log10([438.9 418.9 469.0 512.7]);
E1_8(6:9) = [0.043 0.540 0.258 0.088];


figure(103)
%Subj1 AnthonyM probe on E8
subplot(121)
% errorbar(0:3,M2_6(8:11),E2_6(8:11),'o-m','linewidth',2,'color',[.5 .5 .5]);
errorbar(0:3,M2_8(8:11)-M2_8(8),E2_8(8:11),'o-k','linewidth',2.5); hold on;
% plot(0:4,T2(8:12),'vb',0:4,C2(8:12),'^b','markersize',5); 
% line([0 4],[Lp(6) Lp(6)],'linestyle','--','color',[.5 .5 .5]);
% line([0 4],[Lp(8)-M2_8(8) Lp(8)-M2_8(8)],'linestyle','--','color','k'); 
hold off;
text(4-.1,52.8,'*','color','k','fontsize',18);
% text(-.5,54.5,'E8','edgecolor','k'); %text(-.5,53.8,'E6','edgecolor',[.5 .5 .5]);
set(gca,'XTick',0:4,'xlim',[-1 5],'ylim',[-1 2]); 
% title('AB1 - Adaptively changing MASKER level'); 
ylabel('Masking Re:\Deltax=0 [dB]');  
xlabel({'\Deltax'}); 
%Subj1 AnthonyM single on E8
subplot(122)
% plot(0:3,T1(8:11),'>b',0:3,C1(8:11),'^b',0:3,T1(8:-1:5),'<b',0:3,C1(8:-1:5),'vb','markersize',5); hold on;
line([0 4],[Tp(8) Tp(8)],'linestyle','--','color',[.5 .5 .5]);
line([0 4],[Cp(8) Cp(8)],'linestyle','--','color',[.5 .5 .5]); 
errorbar(0:3,-X2_8(8:11)+X2_8(8),S2_8(8:11),'o-k','linewidth',2.5); hold off;
% text(-.5,54.5,'E8','edgecolor','k');
text(3,Cp(8)-.2-X2_8(8),'Cprobe','color',[.5 .5 .5]);
text(3,Tp(8)+.15-X2_8(8),'Tprobe','color',[.5 .5 .5]);
set(gca,'XTick',0:4,'xlim',[-1 5],'ylim',[-1 2]); 
% title('AB1 - Adaptively changing PROBE level'); ylabel('Amount of masking at probe place [dB]'); 
xlabel({'\Deltax'}); 
% h_legend=legend('T1 basal','C1 basal','T1 apical','C1 apical','location','southeast'); 
% set(h_legend,'FontSize',8,'Position',[0.8679 0.1317 0.1275 0.1783]);
set(gcf,'position',[100 10 800 350]);
%%
plot([331 319 313 358 374 314 411])
 

%% AB SUbj1 (AM) 4 plots
% close all; 
clc; clear all;

Tp([6 8]) = 20*log10([320.0 330.0]);
Cp([6 8]) = 20*log10([273.7 560.0]); 

Lp([6 8]) = 20*log10([345.0 368.0]);

T2(8:12) = 20*log10([240.0 233.7 245.1 261.8 269.4]);
C2(8:12) = 20*log10([390.0 390.0 390.0 405.0 440.0]); %C8,8?

T1(5:11) = 20*log10([269.8 252.9 248.6 294.6 311.5 323.5 351.9]);
C1(5:11) = 20*log10([450.0 450.0 440.0 480.0 490.0 530.0 510.0]); 

M2_8(8:11) = 20*log10([314.7 290.0 313.8 376.0]);
E2_8(8:11) = [0.044 0.192 0.137 0.000];

M2_6(8:11) = 20*log10([290.2 307.4 356.2 379.6]);
E2_6(8:11) = [0.295 0.258 0.407 0.398];

X2_8(8:11) = 20*log10([408.3 400.7 385.4 351.9]);
S2_8(8:11) = [0.254 0.341 0.185 0.133];

M1_8(6:9) = 20*log10([438.9 418.9 469.0 512.7]);
E1_8(6:9) = [0.043 0.540 0.258 0.088];


figure(104)
%Subj1 AnthonyM probe on E8
subplot(121)
% errorbar(0:3,M2_6(8:11),E2_6(8:11),'o-m','linewidth',2,'color',[.5 .5 .5]);
errorbar(0:3,M2_8(8:11),E2_8(8:11),'o-k','linewidth',2.5); hold on;
plot(0:4,T2(8:12),'vb',0:4,C2(8:12),'^b','markersize',5); 
% line([0 4],[Lp(6) Lp(6)],'linestyle','--','color',[.5 .5 .5]);
line([0 4],[Lp(8) Lp(8)],'linestyle','--','color','k'); hold off;
text(4-.1,52.8,'*','color','k','fontsize',18);
% text(-.5,54.5,'E8','edgecolor','k'); %text(-.5,53.8,'E6','edgecolor',[.5 .5 .5]);
set(gca,'XTick',0:4,'xlim',[-1 5],'ylim',[47.0 55.5]); 
% title('AB1 - Adaptively changing MASKER level'); 
ylabel('Masker level at masked Threshold [dB]'); xlabel('\Deltax'); 

%Subj1 AnthonyM single on E8
subplot(122)
plot(0:3,T1(8:11),'>b',0:3,C1(8:11),'^b',0:3,T1(8:-1:5),'<b',0:3,C1(8:-1:5),'vb','markersize',5); hold on;
line([0 4],[Tp(8) Tp(8)],'linestyle','--','color',[.5 .5 .5]);
line([0 4],[Cp(8) Cp(8)],'linestyle','--','color',[.5 .5 .5]); 
errorbar(0:3,X2_8(8:11),S2_8(8:11),'o-k','linewidth',2.5); hold off;
% text(-.5,54.5,'E8','edgecolor','k');
text(3,Cp(8)-.2,'Cprobe','color',[.5 .5 .5]);
text(3,Tp(8)+.15,'Tprobe','color',[.5 .5 .5]);
set(gca,'XTick',0:4,'xlim',[-1 5],'ylim',[47.0 55.5]); 
% title('AB1 - Adaptively changing PROBE level'); 
ylabel('Probe level at masked Threshold [dB]'); xlabel('\Deltax'); 
% h_legend=legend('T1 basal','C1 basal','T1 apical','C1 apical','location','southeast'); 
% set(h_legend,'FontSize',8,'Position',[0.8679 0.1317 0.1275 0.1783]);
set(gcf,'position',[100 10 800 350]);



%% AB SUbj1 (AM) 4 plots

Tp([6 8]) = 20*log10([320.0 330.0]);
Cp([6 8]) = 20*log10([273.7 560.0]); 

Lp([6 8]) = 20*log10([345.0 368.0]);

T2(8:12) = 20*log10([240.0 233.7 245.1 261.8 269.4]);
C2(8:12) = 20*log10([390.0 390.0 390.0 405.0 440.0]); %C8,8?

T1(5:11) = 20*log10([269.8 252.9 248.6 294.6 311.5 323.5 351.9]);
C1(5:11) = 20*log10([450.0 450.0 440.0 480.0 490.0 530.0 510.0]); 

M2_8(8:11) = 20*log10([314.7 290.0 313.8 376.0]);
E2_8(8:11) = [0.044 0.192 0.137 0.000];

M2_6(8:11) = 20*log10([290.2 307.4 356.2 379.6]);
E2_6(8:11) = [0.295 0.258 0.407 0.398];

X2_8(8:11) = 20*log10([408.3 400.7 385.4 351.9]);
S2_8(8:11) = [0.254 0.341 0.185 0.133];

M1_8(6:9) = 20*log10([438.9 418.9 469.0 512.7]);
E1_8(6:9) = [0.043 0.540 0.258 0.088];
figure(105)


%Subj1 AnthonyM dual on E8
subplot(211)
errorbar(0:3,M2_6(8:11),E2_6(8:11),'o-m','linewidth',2,'color',[.5 .5 .5]);hold on;
errorbar(0:3,M2_8(8:11),E2_8(8:11),'o-k','linewidth',2.5); 
plot(0:4,T2(8:12),'vb',0:4,C2(8:12),'^b','markersize',5); 
line([-0.5 4.5],[Lp(6) Lp(6)],'linestyle','--','color',[.5 .5 .5]);
line([-0.5 4.5],[Lp(8) Lp(8)],'linestyle','--','color','k'); hold off;
text(4-.1,52.8,'*','color','k','fontsize',18);
text(-.5,54.5,'E8','edgecolor','k'); text(-.5,53.8,'E6','edgecolor',[.5 .5 .5]);
set(gca,'XTick',0:4,'xlim',[-1 5],'ylim',[47.0 55.5]); title('Dual Masker'); 
ylabel('Masker level at masked Threshold [dB]'); xlabel('\Deltax'); 

%Subj1 AnthonyM single on E8
subplot(212)
plot(-2:1,T1(6:9),'vb',-2:1,C1(6:9),'^b','markersize',5); hold on;
line([-2.5 2.5],[20*log10(350) 20*log10(350)],'linestyle','--','color',[.5 .5 .5]); 
errorbar(-2:1,M1_8(6:9),E1_8(6:9),'o-k','linewidth',2.5); hold off;
text(-2.5,54.5,'E8','edgecolor','k');
set(gca,'XTick',-2:1,'xlim',[-3 3],'ylim',[47.0 55.5]); title('Single Masker'); 
ylabel('Probe level at masked Threshold [dB]'); xlabel('\Deltax'); 
h_legend=legend('T1 basal','C1 basal','location','southeast'); 
% set(h_legend,'FontSize',8,'Position',[0.8679 0.1317 0.1275 0.1783]);
set(gcf,'position',[100 10 350 800]);


%% Table 2 in the MS, values in dB
A=[ 323	368	323	358	408	358
    289	368	289	313	401	 374
    318	368	318	319	385	314
    380	368	380	331	352	411];
round(10*20*log10(A))/10

