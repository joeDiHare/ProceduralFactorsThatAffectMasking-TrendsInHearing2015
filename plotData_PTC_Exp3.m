%%% Exp 3 in PTC article
%% plot IPG investigation Oct.14
clc; clear all; close all; 
facText=[25 50 120 240 380 750 1400 2900];
IPGsText=[30 60 120 240 480 960 1920 7910];


res_all{1}=20*log10(9.4.*[  84.9	84.3	83.1	78.6 %SM
                            77.5	76.7	79.1	76.2
                            70.9	71.4	71.5	74.6
                            60.4	70.4	68.4	70.5
                            68.4	53.9	68.4	67.9
                            66.9	66.4	69.1	66.7
                            73.5	72.5	72.0	71.50
                            67.2	69.7	68.0	64.9
                            67.2	67.6	67.4	68.5
                            66.2	70.0	69.7	69.5
                            67.4	69.7	68.5	69
                            67.2	67.1	66.2	67.8
                            71.2	71.2	75.4	76.4
                            76.7	76.0	82.5	81.9]);
res_all{2}=20*log10(9.4.*[  72.4	74.4	68.6	69.3 %AL
                            69.6	68.9	68.5	70.8
                            68.6	70.1	68.2	69.8
                            65.9	66.4	65.6	66.1
                            64.5	67.2	66.7	66.8
                            62.2	62.2	63.4	65.9
                            64.1	65.3	65.6	65.9
                            63.0	64.3	64.6	63.2
                            63.8	64.1	64.8	64.1
                            63.2	62.5	62.5	63.0
                            62.5	61.5	61.5	62.0
                            64.1	63.4	63.9	63.8
                            74.4	72.6	74.2	73.1
                            77.2	76.5	75.4	75]);
res_all{3}=20*log10(9.4.*[  46.4	48.6	45.6	46.9 %JT
                            47.3	47.0	46.1	44.6
                            44.6	43.0	43.3	43.0
                            40.3	39.7	40.6	40.4
                            39.7	39.5	38.6	38.3
                            37.9	38.4	39.1	38.7
                            42.2	45.3	42.1	42.10
                            42.4	46.5	48.0	46.90
                            43.4	46.2	45.2	48.80
                            46.8	46.4	50.4	50.10
                            38.2	38.4	38.4	38.30
                            38.3	39.6	39.6	38.50
                            50.0	52.0	48.6	48.40
                            53.5	53.8	53.2	53.30]);
res_all{4}=20*log10(9.4.*[  70.0    71.7    67.9    70.6 %CL
                            68.9    68.9    67.2    NaN
                            65.3    64.8    NaN     NaN
                            63.2    63.4    NaN     NaN
                            62.3    61.6    NaN     NaN
                            60.7    60.3    61.6    NaN
                            60.3    61.0    61.4    NaN                    
                            63.2    63.0    NaN     NaN
                            61.6    66.2    62.1    NaN
                            63.9    62.1    NaN     NaN
                            60.7    61.6    NaN     NaN
                            65.3    62.3    NaN     NaN
                            70.9    68.1    NaN     NaN
                            75.4    76.5    74.8    NaN]);
IPGs=[30 60 120 240 480 1920 7910]; IPGp=100*IPGs./8000;
IPGs=[30 60 120 240 480 960 1920 4000 3880 5960 6920 7400 7820 7850]; IPGp=100*IPGs./8000;
IPGs60=[30 60 120 240 480 960 1920  4000 3768 5840 6800 7280 7700 7730]; 


% %% all expressed as amount of masking
% %% plot IPG investigation Oct.14
clc; close(figure(516));%clear all;%  
facText =[25 50 120 240 380 750 1400 2900];
IPGsText=[30 60 120 240 480 960 1920 7910];

figure(516)
h=tight_subplot(2,2,[.00 .00],[.1 .05],[.16 .15]);

% SM
axes(h(1));
res = res_all{1}; nres=nanmean(res,2); nres = nres - nres(8); eres=nanstd(res');
IPGs=[30 60 120 240 480 1920 7910]; IPGp=100*IPGs./8000;
IPGs=[30 60 120 240 480 960 1920 4000 3880 5960 6920 7400 7820 7850]; IPGp=100*IPGs./8000;
probLv=20*log10(9.4*65);
errorbar(IPGs(1:7),nres(1:7)',eres(1:7)','o-k','linewidth',2,'markerfacecolor','k');     hold on;
errorbar(8000-IPGs(end:-1:10)-120,nres(end:-1:10)',eres(end:-1:10)','o-r','linewidth',2,'markerfacecolor','w');
errorbar(4e3,nres(8)',eres(8)','o-k','linewidth',2,'markerfacecolor','k');     
errorbar(4e3,nres(9)',eres(9)','o-r','linewidth',2,'markerfacecolor','w');
line([30 4000],[0 0],'color',[.5 .5 .5],'linestyle','--'); hold off;
set(gca,'XTick',IPGs(1:8),'ytick',-1.5:.5:2,'xticklabel',''); 
text(1e3,2,'ME1','fontweight','bold');%title('ME3 (SM); single masker, interleaved masking, 125pps'); 
for n=1, text(facText(n),54.7,['-',num2str(IPGsText(n))],'fontsize',10,'color','r'); end; hold off; %[1:5 7]
ylabel('MLTs  {\itre:} MPG=4000\mus [dB]'); xlabel('MPG [\mus]'); 
axis([15 8000 -2.0 2.5]); set(gca, 'XScale', 'log');
% AL
axes(h(2));
res = res_all{2}; nres=nanmean(res,2); nres = nres - nres(8); eres=nanstd(res');
IPGs=[30 60 120 240 480 960 1920  4000 3768 5840 6800 7280 7700 7730]; IPGp=100*IPGs./8000;
probLv=20*log10(9.4*65);
errorbar(IPGs(1:7),nres(1:7)',eres(1:7)','o-k','linewidth',2,'markerfacecolor','k');     hold on;
errorbar(8000-IPGs(end:-1:10)-240,nres(end:-1:10)',eres(end:-1:10)','o-r','linewidth',2,'markerfacecolor','w');
errorbar(4e3,nres(8)',eres(8)','o-k','linewidth',2,'markerfacecolor','k');     
errorbar(4e3,nres(9)',eres(9)','o-r','linewidth',2,'markerfacecolor','w');
line([30 4000],[0 0],'color',[.5 .5 .5],'linestyle','--')
set(gca,'XTick',IPGs(1:8),'ytick',-1.5:.5:2,'yticklabel','','xticklabel',''); 
text(1e3,2,'ME2','fontweight','bold');%title('ME4 (AL); single masker, interleaved masking, 125pps'); 
for n=[1:2 5:8], text(facText(n),55.0,['-',num2str(IPGs(n))],'fontsize',10,'color','r'); end
xlabel('MPG [\mus]');
axis([15 8000 -2.0 2.5]); set(gca, 'XScale', 'log');
%%% JT
axes(h(3));
res = res_all{3}; nres=nanmean(res,2); nres = nres - nres(8); eres=nanstd(res');
IPGs=[30 60 120 240 480 960 1920 4000 3880 5960 6920 7400 7820 7850]; IPGp=100*IPGs./8000;
probLv=20*log10(9.4*39);
errorbar(IPGs(1:7),nres(1:7)',eres(1:7)','o-k','linewidth',2,'markerfacecolor','k');     hold on;
errorbar(8000-IPGs(end:-1:10)-120,nres(end:-1:10)',eres(end:-1:10)','o-r','linewidth',2,'markerfacecolor','w');
errorbar(4e3,nres(8)',eres(8)','o-k','linewidth',2,'markerfacecolor','k');     
errorbar(4e3,nres(9)',eres(9)','o-r','linewidth',2,'markerfacecolor','w');
line([30 4000],[0 0],'color',[.5 .5 .5],'linestyle','--'); hold off;
set(gca,'XTick',IPGs(1:8),'ytick',-1.5:.5:2); 
text(1e3,2,'ME3','fontweight','bold');%title('ME5 (JT); single masker, interleaved masking, 125pps'); 
for n=[1 2 5:8], text(facText(n),50.7,['-',num2str(IPGsText(n))],'fontsize',10,'color','r'); end
ylabel('MLTs  {\itre:} MPG=4000\mus [dB]'); xlabel('MPG [\mus]'); 
axis([15 8000 -2.0 2.5]); set(gca, 'XScale', 'log');
%%% CL
axes(h(4));
res = res_all{4}; nres=nanmean(res,2); nres = nres - nres(8); eres=nanstd(res');
IPGs=[30 60 120 240 480 960 1920 4000 3880 5960 6920 7400 7820 7850]; IPGp=100*IPGs./8000;
probLv=20*log10(9.4*60);
errorbar(IPGs(1:7),nres(1:7)',eres(1:7)','o-k','linewidth',2,'markerfacecolor','k');     hold on;
errorbar(8000-IPGs(end:-1:10)-120,nres(end:-1:10)',eres(end:-1:10)','o-r','linewidth',2,'markerfacecolor','w');
errorbar(4e3,nres(8)',eres(8)','o-k','linewidth',2,'markerfacecolor','k');     
errorbar(4e3,nres(9)',eres(9)','o-r','linewidth',2,'markerfacecolor','w');
line([30 4000],[0 0],'color',[.5 .5 .5],'linestyle','--'); hold off;
set(gca,'XTick',IPGs(1:8),'ytick',-1.5:.5:2,'yticklabel',''); 
text(1e3,2,'ME5','fontweight','bold');%title('ME5 (JT); single masker, interleaved masking, 125pps'); 
for n=[1 5:7], text(IPGsText(n)-facText(n),50.7,['-',num2str(IPGsText(n))],'fontsize',10,'color','r'); end
for n=[1:2 5:8], text(facText(n),54.35,['-',num2str(IPGs(n))],'fontsize',10,'color','r'); end
 xlabel('MPG [\mus]');
axis([15 8000 -2.0 2.5]); set(gca, 'XScale', 'log');
set(gcf,'position',[100 100 800 550]);
print(gcf, '-dpng', '-r200', 'Fig_7_PTC.png');


%%
figure(517)
h=tight_subplot(2,2,[.01 .06],[.1 .05],[.16 .15]);


%%% S1
axes(h(1));
res = res_all{1}; nres=nanmean(res,2); eres=nanstd(res');
probLv=20*log10(9.4*65);
errorbar(IPGs(1:7),nres(1:7)',eres(1:7)','o-k','linewidth',2);     hold on;
errorbar(8000-IPGs(end:-1:10)-120,nres(end:-1:10)',eres(end:-1:10)','o-r','linewidth',2);
errorbar(4e3,nres(8)',eres(8)','o-k','linewidth',2);     
errorbar(4e3,nres(9)',eres(9)','o-r','linewidth',2);
line([30 8000],[probLv probLv],'color',[.5 .5 .5],'linestyle','--'); hold off;
set(gca,'XTick',IPGs(1:8)); 
text(4000,58.0,'S3');%title('ME3 (SM); single masker, interleaved masking, 125pps'); 
for n=[1 2 5:8], text(facText(n),54.7,['-',num2str(IPGsText(n))],'fontsize',10,'color','r'); end
ylabel({'Masker level at masked threshold' '( \Deltax = 0 ) [dB]'}); xlabel('masker-probe gap ( \mus )'); 
axis([10 8000 54.5 58.5]); set(gca, 'XScale', 'log');

% AL
axes(h(2));
res = res_all{2}; nres=nanmean(res,2); eres=nanstd(res');
probLv=20*log10(9.4*65);
errorbar(IPGs60(1:7),nres(1:7)',eres(1:7)','o-k','linewidth',2);     hold on;
errorbar(8000-IPGs60(end:-1:10)-240,nres(end:-1:10)',eres(end:-1:10)','o-r','linewidth',2);
errorbar(4e3,nres(8)',eres(8)','o-k','linewidth',2);     
errorbar(4e3,nres(9)',eres(9)','o-r','linewidth',2);
line([30 4000],[probLv probLv],'color',[.5 .5 .5],'linestyle','--')
set(gca,'XTick',IPGs60(1:8)); 
text(4000,57.1,'S4');%title('ME4 (AL); single masker, interleaved masking, 125pps'); 
for n=[1:2 5:8], text(facText(n),55.0,['-',num2str(IPGs60(n))],'fontsize',10,'color','r'); end
xlabel('masker-probe gap ( \mus )'); 
axis([10 8000 54.9 57.5]); set(gca, 'XScale', 'log');

%%% JT
axes(h(3));
res = res_all{3}; nres=nanmean(res,2); eres=nanstd(res');
probLv=20*log10(9.4*39);
errorbar(IPGs(1:7),nres(1:7)',eres(1:7)','o-k','linewidth',2);     hold on;
errorbar(8000-IPGs(end:-1:10)-120,nres(end:-1:10)',eres(end:-1:10)','o-r','linewidth',2);
errorbar(4e3,nres(8)',eres(8)','o-k','linewidth',2);     
errorbar(4e3,nres(9)',eres(9)','o-r','linewidth',2);
line([30 8000],[probLv probLv],'color',[.5 .5 .5],'linestyle','--'); hold off;
set(gca,'XTick',IPGs(1:8)); 
text(4000,54,'S5');%title('ME5 (JT); single masker, interleaved masking, 125pps'); 
for n=[1 2 5:8], text(facText(n),50.7,['-',num2str(IPGsText(n))],'fontsize',10,'color','r'); end
ylabel({'Masker level at masked threshold' '( \Deltax = 0 ) [dB]'});xlabel('masker-probe gap ( \mus )'); 
axis([10 8000 50.5 54.5]); set(gca, 'XScale', 'log');

%%% CL
axes(h(4));
res = res_all{4}; nres=nanmean(res,2); eres=nanstd(res');
probLv=20*log10(9.4*60);
errorbar(IPGs(1:7),nres(1:7)',eres(1:7)','o-k','linewidth',2);     hold on;
errorbar(8000-IPGs(end:-1:10)-120,nres(end:-1:10)',eres(end:-1:10)','o-r','linewidth',2);
errorbar(4e3,nres(8)',eres(8)','o-k','linewidth',2);     
errorbar(4e3,nres(9)',eres(9)','o-r','linewidth',2);
line([30 8000],[probLv probLv],'color',[.5 .5 .5],'linestyle','--'); hold off;
set(gca,'XTick',IPGs(1:8)); 
text(4000,57,'S9');%title('ME5 (JT); single masker, interleaved masking, 125pps'); 
for n=[1 5:7], text(IPGsText(n)-facText(n),50.7,['-',num2str(IPGsText(n))],'fontsize',10,'color','r'); end
for n=[1:2 5:8], text(facText(n),54.35,['-',num2str(IPGs(n))],'fontsize',10,'color','r'); end
xlabel('masker-probe gap ( \mus )'); 
axis([10 8000 54.2 57.5]); set(gca, 'XScale', 'log');


% %%% modelling
% subplot(224)
% t=linspace(31,8000,5e3); for i=1:length(IPGs), ind(i)=argmin(abs(t-IPGs(i))); end
% nf = 10*(30./(t)).^1; 
% s=100; mu=30; nf = 1000*exp(-.5 * ((t - mu)./s) .^ 2) ./ (s * sqrt(2*pi)); 
% rf = -1./(1+10.^(0.05*(t-1200)));
% rb = -.5./(1+10.^(0.01*(t-1000)));
% plot(t,nf,'m--',t,rf,'b-.','linewidth',2);          hold on;
% comb_f=nf+rf; comb_b=nf+rb; 
% plot(IPGs,comb_f(ind),'o-','linewidth',3,'color','k'); 
% % plot(IPGs,comb_b(ind),'x-','linewidth',3,'color',[.5 .5 .5]); 
% % plot(t,rb,'k.','markersize',4);          
% plot([1 8000],[0 0],'--','color',[.7 .7 .7]);   hold off;
% IPGs=[30 60 120 240 480 960 1920];
% set(gca,'XTick',IPGs); 
% set(gca, 'XScale', 'log');axis([10 8000 -1.6 5.6]);
% xlabel 'masker-probe gap /mus'; ylabel 'gain on probe detection'; title 'simple model';
% legend('excitation (facilitation)','inhibition (refractoriness)','effects combined');

set(gcf,'position',[100 100 1050 600]);

%%
clc; %clear all; 
facText=[25 50 120 240 380 750 1400 2900];
IPGsText=[30 60 120 240 480 960 1920 7910];

% SM
res = res_all{1}; nres1=nanmean(res,2); eres1=nanstd(res'); probLv1=20*log10(9.4*65);
% AL
res = res_all{2}; nres2=nanmean(res,2); eres2=nanstd(res'); probLv3=20*log10(9.4*65);
% JT
res = res_all{3}; nres3=nanmean(res,2); eres3=nanstd(res'); probLv2=20*log10(9.4*39);
% CL
res = res_all{4}; nres4=nanmean(res,2); eres4=nanstd(res'); probLv4=20*log10(9.4*60);
                
% nres = nanmean([nres1 nres2 nres3 nres4],2);
nres = nanmean([nres1-nres1(8) nres2-nres2(8) nres3-nres3(8) nres4-nres4(8)],2);%
eres = nanstd([nres1-nres1(8) nres2-nres2(8) nres3-nres3(8) nres4-nres4(8)],[],2)./4;%
probLv = mean([probLv1 probLv2 probLv3 probLv4],2);

figure(181)
subplot(2,6,8:11);%subplot(2,6,2:5)
IPGs=[30 60 120 240 480 960 1920 4000 3880 5960 6920 7400 7820 7850]; IPGp=100*IPGs./8000;
errorbar(IPGs(1:8),nres(1:8)',eres(1:8)','o-k','linewidth',2,'markerfacecolor','k');     hold on;
errorbar(8000-IPGs(end:-1:9)-120,nres(end:-1:9)',eres(end:-1:9)','o-r','linewidth',2,'markerfacecolor','r');
line([30 8000],[probLv probLv],'color',[.5 .5 .5],'linestyle','--'); 
set(gca,'XTick',IPGs(1:8)); 
text(4000,57,'S9');%title('ME5 (JT); single masker, interleaved masking, 125pps'); 
for n=[1 5:7], text(IPGsText(n)-facText(n),50.7,['-',num2str(IPGsText(n))],'fontsize',10,'color','r'); end
for n=[1:2 5:8], text(facText(n),53.6,['-',num2str(IPGs(n))],'fontsize',10,'color','r'); end
xlabel('masker-probe gap ( \mus )'); ylabel({'Relesae of masking' 'MLTs  re:MLT  at  4ms (dB)'}); title 'Average Normalised Data (N=4)';
plot([1 8000],[0 0],'--','color',[.7 .7 .7]);   hold off;
axis([10 8000 -1 2]); set(gca, 'XScale', 'log');


% %% % modelling - forward IPGs
subplot(2,6,1:3);%subplot(2,6,7:9)
% fun=@(A,X)(A(1)./(X.^A(2)) + A(3)./( 1+(10.^(A(4)*(X-A(5)))))  ); A0 = [30 1 -1 .05 1200 -1]; % Initial values
fun=@(A,X)(A(1)*exp(-.5 * ((X - A(2))./A(3)) .^ 2) ./ (A(3) * sqrt(2*pi)) ...
         + A(4)./( 1+(10.^(A(5)*(X-A(6)))))  ); A0 = [192610   -1600    529   -0.6013    0.0179    1913]; % Initial values
x = IPGs(1:8); y = nres(1:8)'; ey = eres(1:8);
A_fit = real(nlinfit(x, y, fun, A0)); A_fit=[980000 -1650 473 -0.5 .0025 1355];
x_range=linspace(min(IPGs),max(IPGs),1e3); y_sigmo=real(fun(A_fit,x_range));
nf=real(fun([A_fit(1) A_fit(2) A_fit(3) 0 1 1],x_range));
rf=real(fun([0 1 1 A_fit(4) A_fit(5) A_fit(6)],x_range));
% [param,h,x_range,y_sigmo]=sigm_fit(PC_durats(1:end,1),PC_durats(1:end,2),[],[.1 .1 .1 .1 .1 .9],0);
for i=1:length(x), ind(i)=argmin(abs(x_range-x(i))); end
plot(x,y_sigmo(ind),'-','color',[.5 .5 .5],'linewidth',5); hold on;
errorbar(x,y,ey,'o-k','linewidth',2,'markerfacecolor','k');
plot(x_range(1:509),nf(1:509),'--b',x_range(1:509),rf(1:509),'--g','linewidth',2);        
plot([1 8000],[0 0],'--','color',[.7 .7 .7]);   hold off;
set(gca,'XTick',x, 'XScale', 'log'); axis([10 8000 -1 2]);
xlabel 'masker-probe gap ( \mus )'; ylabel({'Predicted MLT re: 4 ms (dB)'}); title 'Data Modelling - forward IPGs';
legend('model','data','excitation (facilitation)','inhibition (refractoriness)');
rmse(y,y_sigmo(ind))  % err:0.0532


%% % modelling - backward IPGs
subplot(2,6,4:6)% subplot(2,6,10:12)
% fun=@(A,X)(A(1)*exp(-.5 * ((X - A(2))./A(3)) .^ 2) ./ (A(3) * sqrt(2*pi)) ...
%          + A(4)./( 1+(10.^(A(5)*(X-A(6)))))  ); A0 = [192610   -1600    529   -0.465    0.0179    1913]; % Initial values
x = 8000-IPGs(end:-1:9)-120; y = nres(end:-1:9)'; ey = eres(end:-1:9)';
% fun=@(A,X)(A(1)*exp(-.5 * ((X - A(2))./A(3)) .^ 2) ./ (A(3) * sqrt(2*pi)) ...
%          + A(4)./( 1+(10.^(A(5)*(X-A(6)))))  ); A0 = [192610   -1600    529   -0.6013    0.0179    1913]; % Initial values
% A_fit = real(nlinfit(x, y, fun, A0)); A_fit=[980000 -1650 473 -0.2 .0025 1355];
fun=@(A,X)(A(1)*exp(-.5 * ((X - A(2))./A(3)) .^ 2) ./ (A(3) * sqrt(2*pi)) ...
         + A(4)./( 1+(10.^(A(5)*(X-A(6)))))  ); A0 = [192610   -1600    529   -0.6013    0.0179    1913]; % Initial values
A_fit = real(nlinfit(x, y, fun, A0)); A_fit=[980000 -1650 473 -0.5 .0025 1355];

x_range=linspace(min(IPGs),max(IPGs),1e3); y_sigmo=real(fun(A_fit,x_range));
nf=real(fun([A_fit(1) A_fit(2) A_fit(3) 0 1 1],x_range));
rf=real(fun([0 1 1 A_fit(4) A_fit(5) A_fit(6)],x_range));
ind=[]; for i=1:length(x), ind(i)=argmin(abs(x_range-x(i))); end
plot(x,y_sigmo(ind),'-','color',[.9 .5 .5],'linewidth',5); hold on;
errorbar(x,y,ey,'o-r','linewidth',2,'markerfacecolor','r');
plot(x_range(1:509),nf(1:509),'--b',x_range(1:509),rf(1:509),'--m','linewidth',2);        
plot([1 8000],[0 0],'--','color',[.7 .7 .7]);   hold off;
set(gca,'XTick',x, 'XScale', 'log'); axis([10 8000 -1 2]);
title 'Data Modelling - backward IPGs'; xlabel 'masker-probe gap ( \mus )'; %ylabel({'Predicted MLT re: 4 ms (dB)'}); 
legend('model','data','excitation (facilitation)','inhibition (refractoriness)');
rmse(y,y_sigmo(ind))
%% % modelling - backward IPGs CHANGING REF ONLY
subplot(2,6,4:6)% subplot(2,6,10:12)
% fun=@(A,X)(A(1)*exp(-.5 * ((X - A(2))./A(3)) .^ 2) ./ (A(3) * sqrt(2*pi)) ...
%          + A(4)./( 1+(10.^(A(5)*(X-A(6)))))  ); A0 = [192610   -1600    529   -0.465    0.0179    1913]; % Initial values
x = 8000-IPGs(end:-1:9)-120; y = nres(end:-1:9)'; ey = eres(end:-1:9)';
fun=@(A,X)(980000*exp(-.5 * ((X +1650)./473) .^ 2) ./ (473 * sqrt(2*pi)) ...
         + A(1)*exp(-.5 * ((X - A(2))./A(3)) .^ 2) ./ (A(3) * sqrt(2*pi))  ); 
     A0 = [-0.6013    0.0179    1913]; % Initial values
A_fit = real(nlinfit(x, y, fun, A0)); %A_fit=[980000 -1650 473 -0.2 .0025 1355];
A_fit=[-550 700 450];
x_range=linspace(min(IPGs),max(IPGs),1e3); y_sigmo=real(fun(A_fit,x_range));
% nf=real(fun([A_fit(1) A_fit(2) A_fit(3) 0 1 1],x_range));
fun_fr=@(A,X)(A(1)*exp(-.5 * ((X - A(2))./A(3)) .^ 2) ./ (A(3) * sqrt(2*pi))  ); 
rf=real(fun_fr(A_fit,x_range));
ind=[]; for i=1:length(x), ind(i)=argmin(abs(x_range-x(i))); end
plot(x,y_sigmo(ind),'-','color',[.9 .5 .5],'linewidth',5); hold on;
errorbar(x,y,ey,'o-r','linewidth',2,'markerfacecolor','r');
plot(x_range(1:509),nf(1:509),'--b',x_range(1:509),rf(1:509),'--m','linewidth',2);        
plot([1 8000],[0 0],'--','color',[.7 .7 .7]);   hold off;
set(gca,'XTick',x, 'XScale', 'log'); axis([10 8000 -1 2]);
title 'Data Modelling - backward IPGs'; xlabel 'masker-probe gap ( \mus )'; %ylabel({'Predicted MLT re: 4 ms (dB)'}); 
legend('model','data','excitation (facilitation)','inhibition (refractoriness)');
rmse(y,y_sigmo(ind))

set(gcf,'position',[100 100 950 600]);


