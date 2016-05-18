%%% Exp. 4 in PTC article
%% confusion effect PSY Nov14
figure(431); 
clear all; %clc;
interc=51/100;
semilogx(0,0,'dk',0,0,'ok',0,0,'sk',0,0,'^k','linewidth',2,'markerfacecolor',grey); hold on;%,0,0,'.r'
semilogx(0,0,'dk',0,0,'ok',0,0,'sk',0,0,'^k','linewidth',2,'markerfacecolor','m');
legend('ME1','ME2','ME3','ME5','location','northwest');%,...
%        'rem. mask ME1','rem. mask ME2','rem. mask ME3','rem. mask ME5','location','northwest');%,'remote masker'

%%% data
PC_durats_all{1}=[0.32      0.99     %ME3 [SM]
             0.1600    0.935
             0.0800    0.665
             0.0400    0.555
             0.02      0.565
             0.01      0.507
             0.005     0.515];
PC_durats_all{2}=[0.32      0.92 %AL
             0.1600    0.7900
             0.0800    0.6500
             0.0400    mean([0.5000 .5])
             0.02      mean([0.44 0.52])
             0.01      .5
             0.005     .5];
PC_durats_all{3}=[0.32      1 %JT
             0.1600    .98
             0.0800    1
             0.0400    mean([0.66 .72])
             0.03      0.62
             0.02      .5
             0.005     .5]; 
PC_durats_all{4}=[0.32      1 %CL
             0.1600    1
             0.0800    1
             0.0400    mean([0.82 .78])
             0.02      .66
             0.01      mean([.52 .54])
             0.005     mean([.52 .48])]; 

%%%ME3 [SM]
PC_durats = PC_durats_all{1};
[x,b]=sort(PC_durats(:,1)); y=PC_durats(b,2);
%%% estimage sigmoidal fit
sigfunc=@(A,X)(.5 + (max(y)-.5)./(1+exp(-A(1)*(X-A(2))))); A0 = [50 .5]; % Initial values
sigfunc=@(A,X)(.5 + (max(y)-.5)./( 1+(X./A(1)).^A(2) ) );  A0 = [.05 1]; % Initial values
A_fit = real(nlinfit(x, y, sigfunc, A0)); %A_fit=[-0.116   -2.38896];
x_range=linspace(min(x),max(x),1e3); y_sigmo=normaNr(real(sigfunc(A_fit,x_range)),.5,max(y));
% [param,h,x_range,y_sigmo]=sigm_fit(PC_durats(1:end,1),PC_durats(1:end,2),[],[.1 .1 .1 .1 .1 .9],0);
Y_sig=[x_range(:) y_sigmo(:)];
x_interc=argmin(abs(interc-Y_sig(:,2))); y_interc=(x_range(x_interc));
x_interc75=argmin(abs(.75-Y_sig(:,2)));  y_interc75(1)=(x_range(x_interc75));
for i=1:length(x), Yx(i) = argmin(abs(x_range-x(i))); end; e1=rmse(y',y_sigmo(Yx));
semilogx(x,100*y,'dk',x_range,100*y_sigmo,'k','linewidth',2,'markerfacecolor',grey,'markersize',10);
semilogx([0.02 0.08],[57.5 69.5],'dm',[0.08 0.08],[66.5 69.5],'-.m','linewidth',2,'markersize',10);

%%%ME4 [AL]
PC_durats = PC_durats_all{2};
[x,b]=sort(PC_durats(:,1)); y=PC_durats(b,2);
%%% estimage sigmoidal fit
sigfunc=@(A,X)(.5 + (max(y)-.5)./(1+exp(-A(1)*(X-A(2))))); A0 = [50 .5]; % Initial values
sigfunc=@(A,X)(.5 + (max(y)-.5)./( 1+(X./A(1)).^A(2) ) ); A0 = [.05 1]; % Initial values
A_fit = real(nlinfit(x, y, sigfunc, A0)); A_fit=[-0.116   -2.38896];
x_range=linspace(min(x),max(x),1e3); y_sigmo=normaNr(real(sigfunc(A_fit,x_range)),.5,max(y));
% [param,h,x_range,y_sigmo]=sigm_fit(PC_durats(1:end,1),PC_durats(1:end,2),[],[.1 .1 .1 .1 .1 .9],0);
Y_sig=[x_range(:) y_sigmo(:)];
x_interc=argmin(abs(interc-Y_sig(:,2))); y_interc=(x_range(x_interc));
x_interc75=argmin(abs(.75-Y_sig(:,2)));  y_interc75(2)=(x_range(x_interc75));
for i=1:length(x), Yx(i) = argmin(abs(x_range-x(i))); end; e2=rmse(y',y_sigmo(Yx));
semilogx(x,100*y,'ok',x_range,100*y_sigmo,'k','linewidth',2,'markerfacecolor',grey,'markersize',10); 
semilogx(0.04,100*mean([.5 .54]),'om',[0.04 0.04],100*[y(2) mean([.54 .5])],'-.m','linewidth',2,'markersize',10);
% semilogx(0.08,100*0.735,'or',[0.08 0.08],100*[.735 .78],'--r','markerfacecolor','r');

%%%ME5 [JT]
PC_durats = PC_durats_all{3};
[x,b]=sort(PC_durats(:,1)); y=PC_durats(b,2);
%%% estimage sigmoidal fit
% sigfunc=@(A,X)(0 + (max(y))./(1+exp(-A(1)*(X-A(2))))); A0 = [50 .5]; % Initial values
sigfunc=@(A,X)( (max(y)-.5)./( 1+(X./A(1)).^A(2) ) +.5 ); A0 = [50 1]; % Initial values
A_fit = real(nlinfit(x, y, sigfunc, A0)); A_fit=[0.042   -4.831];
x_range=linspace(min(x),max(x),1e3);  y_sigmo=normaNr(real(sigfunc(A_fit,x_range)),.5,max(y));
[param,h,x_range,y_sigmo]=sigm_fit(PC_durats(1:end,1),PC_durats(1:end,2),[],[.1 .1 .1 .1 .1 .9],0);
y_sigmo=normaNr(y_sigmo,.5,max(y));
Y_sig=[x_range(:) y_sigmo(:)];
x_interc=argmin(abs(interc-Y_sig(:,2))); y_interc=(x_range(x_interc));
x_interc75=argmin(abs(.75-Y_sig(:,2)));  y_interc75(3)=(x_range(x_interc75));
for i=1:length(x), Yx(i) = argmin(abs(x_range-x(i))); end; e3=rmse(y',y_sigmo(Yx));
semilogx(x,100*y,'sk',x_range,100*y_sigmo,'k','linewidth',2,'markerfacecolor',grey,'markersize',10);
semilogx(0.03,100*.57,'sm',[0.03 0.03],100*[y(3) .57],'-.m','linewidth',2,'markersize',10);

%%%ME9 [CL]
PC_durats = PC_durats_all{4};
[x,b]=sort(PC_durats(:,1)); y=PC_durats(b,2);
%%% estimage sigmoidal fit
sigfunc=@(A,X)(0 + (max(y))./(1+exp(-A(1)*(X-A(2))))); A0 = [50 .5]; % Initial values
sigfunc=@(A,X)( (max(y)-.5)./( 1+(X./A(1)).^A(2) ) +.5 ); A0 = [50 1]; % Initial values
A_fit = real(nlinfit(x, y, sigfunc, A0)); A_fit=[0.0299   -3.031];
x_range=linspace(min(x),max(x),1e3);  y_sigmo=normaNr(real(sigfunc(A_fit,x_range)),.5,max(y));
% [param,h,x_range,y_sigmo]=sigm_fit(PC_durats(1:end,1),PC_durats(1:end,2),[],[.1 .1 .1 .1 .1 .9],0);y_sigmo=normaNr(y_sigmo,.5,max(y));
Y_sig=[x_range(:) y_sigmo(:)];
x_interc=argmin(abs(interc-Y_sig(:,2))); y_interc=(x_range(x_interc));
x_interc75=argmin(abs(.75-Y_sig(:,2)));  y_interc75(4)=(x_range(x_interc75));
for i=1:length(x), Yx(i) = argmin(abs(x_range-x(i))); end; e4=rmse(y',y_sigmo(Yx));
semilogx(x,100*y,'^k',x_range,100*y_sigmo,'k','linewidth',2,'markerfacecolor',grey,'markersize',10); 
semilogx(0.02,100*mean([.54 .58]),'^m',[0.02 0.02],100*[y(3) mean([.54 .58])],'-.m','linewidth',2,'markersize',10); hold off;

hold off;
[100.*[e1 e2 e3 e4]; round(1e3.*y_interc75)]
line([1 800].*1e-3, 100.*[.5 .5],'linestyle','--','color',[.5 .5 .5])
xlabel('probe duration [ms]'); ylabel 'percent correct';
% title(['intercept at ',num2str(interc,2),'% is: ',num2str(1000*y_interc),'ms']);
axis([0.003 0.4 100*.42 100*1.05]);
set(gca,'xtick',sort(x),'xticklabels',{'5','10','20','40','80','160','320'});
print(gcf, '-dpng', '-r200', 'Fig_8_PTC.png');

%% last plot Feb.15 on AL
%%%ME4 [AL]
figure(432)
subplot(211)
errorbar(0,62.5,8.2,'ok','markerfacecolor','k','markersize',8,'linewidth',1.5); hold on;
errorbar(1,75.0,5.0,'ok','markerfacecolor','k','markersize',8,'linewidth',1.5);
errorbar(2,73.5,6.9,'ok','markerfacecolor','r','markersize',8,'linewidth',1.5);
errorbar(3,54.5,3.1,'ok','markerfacecolor','k','markersize',8,'linewidth',1.5);
errorbar(4,54.5,2.0,'ok','markerfacecolor','r','markersize',8,'linewidth',1.5); 
plot([-.2 5.2],[50 50],'--','color',grey); hold off;
cond = {'         P                   M             no-Rem    ' ...
        '      P+1dB          M+1dB          no-Rem      ' ...
        '      P+1dB          M+1dB          Rem 11      ' ...
        '          P              M+1dB          no-Rem      ' ... 
        '          P              M+1dB          Rem 11      '};
axis([-.5 4.5 40 85]); set(gca,'xtick',0:4,'xticklabel',cond); fix_xticklabels();
title('Subj.4 [AL] on confusion effect'); ylabel('percent correct'); xlabel('conditions'); 

subplot(212)
errorbar(0,65.0, 2.0,'ok','markerfacecolor','k','markersize',8,'linewidth',1.5); hold on;
errorbar(1,57.0, 8.7,'ok','markerfacecolor','r','markersize',8,'linewidth',1.5);
errorbar(2,61.3,10.9,'ok','markerfacecolor','k','markersize',8,'linewidth',1.5);
errorbar(3,61.0, 9.9,'ok','markerfacecolor','k','markersize',8,'linewidth',1.5);
errorbar(4,50.0, 8.8,'ok','markerfacecolor','k','markersize',8,'linewidth',1.5); 
plot([-.2 5.2],[50 50],'--','color',grey); hold off;
cond = {'         P                  M             no-Rem    ' ...
        '         P                  M             Rem 11      ' ...
        '      P+1dB          M+1dB          no-Rem      ' ...
        '        P             M+0.5dB          no-Rem      ' ... 
        '        P              M+1dB           no-Rem      '};
axis([-.5 4.5 40 85]); set(gca,'xtick',0:4,'xticklabel',cond); fix_xticklabels();
title('Subj.5 [JT] on confusion effect'); ylabel('percent correct'); xlabel('conditions'); 

set(gcf,'position',[70 10 600 900])
