%% validation
close all
clear all

% select reference level and time
z0=200;
t0=[1980 2000];

% coordinates
x=-179.5:179.5;
y=-89.5:89.5;
xm=1:2:359;
ym=-89:2:89;
load z33.mat;
z33=z;
time=[1850:2014]+0.5;
It=find(time>=t0(1)&time<=t0(2));

% open files 
dirs={'/Users/takamitsu3/Dropbox (GaTech)/Research/MI-CMIP6/WOA2018/' ... 
    '/Users/takamitsu3/Dropbox (GaTech)/Research/MI-CMIP6/CMIP6/'};

% load obs
z=ncread([dirs{1},'woa18_all_o00_01.nc'],'depth');
k=find(z==z0);
o2obs=ncread([dirs{1},'woa18_all_o00_01.nc'],'o_an',[1 1 k 1],[360 180 1 1]);
Tobs=ncread([dirs{1},'woa18_decav_t00_01.nc'],'t_an',[1 1 k 1],[360 180 1 1]);
Sobs=ncread([dirs{1},'woa18_decav_s00_01.nc'],'s_an',[1 1 k 1],[360 180 1 1]);

for m=1:4
    model={'CanESM5' 'CNRM-ESM2-1' 'GFDL-CM4' 'IPSL-CM6A-LR'};
    run={'r1i1p1f1' 'r1i1p1f2' 'r1i1p1f1' 'r1i1p1f1'};
    % load model
    k=find(z33==z0);
    tmp=load([dirs{2},'/',model{m},'/o2_',model{m},'_historical_',run{m},'_2x2L33_ann_dft.mat']); 
    tmp2=mean(tmp.o2(:,:,1:k,It),4);
    for i=2:179
        for j=1:90
            if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
            end
        end
    end
    o2mod(:,:,m)=tmp2(:,:,end)*1e3;
    o2mod3d(:,:,:,m)=tmp2*1e3;
    
    clear tmp;
    tmp=load([dirs{2},'/',model{m},'/thetao_',model{m},'_historical_',run{m},'_2x2L33_ann_dft.mat']); 
    tmp2=mean(tmp.thetao(:,:,1:k,It),4);
    for i=2:179
        for j=1:90
            if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
            end
        end
    end
    Tmod(:,:,m)=tmp2(:,:,end);
    Tmod3d(:,:,:,m)=tmp2;
    
    clear tmp;
    tmp=load([dirs{2},'/',model{m},'/so_',model{m},'_historical_',run{m},'_2x2L33_ann_dft.mat']); 
    tmp2=mean(tmp.so(:,:,1:k,It),4);
    for i=2:179
        for j=1:90
            if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
            end
        end
    end
    Smod3d(:,:,:,m)=tmp2*1e3;
    Smod(:,:,m)=tmp2(:,:,end);

end

x(end+1)=180.5;
xm(end+1)=361;
[xo,yo]=meshgrid(x,y);
[xx,yy]=meshgrid(xm,ym);
o2obs(end+1,:)=o2obs(1,:);
Tobs(end+1,:)=Tobs(1,:);
Sobs(end+1,:)=Sobs(1,:);
o2mod(end+1,:,:)=o2mod(1,:,:);
Tmod(end+1,:,:)=Tmod(1,:,:);
Smod(end+1,:,:)=Smod(1,:,:);

addpath ~/MATLAB/m_map
% O2 comparison
figure(1);
subplot1(3,2,'Gap',[0.02 0.02]);
subplot1(6);
m_proj('robinson','clon',-150);
m_pcolor(xo,yo,o2obs');
hold on;
m_pcolor(xo-360,yo,o2obs');
hold off;
shading flat;
caxis([0 350]);
cmp=colormap('jet');
m_grid('xticklabels','','yticklabels','');
m_coast;

for m=1:4
    subplot1(m);
    m_proj('robinson','clon',-150);
    m_pcolor(xx,yy,o2mod(:,:,m)');
    hold on;
    m_pcolor(xx-360,yy,o2mod(:,:,m)');
    hold off;
    shading flat;
    colormap('jet');
    m_grid('xticklabels','','yticklabels','');
    m_coast;
    caxis([0 350]);
end
subplot1(5);
set(gca,'visible','off');
colorbartype([.1 .2 .38 .02],0:10:350,36,[0 350],cmp,0);
set(gca,'xtick',1:5:36);
set(gca,'xticklabel',{'0' '50' '100' '150' '200' '250' '300' '350'});
title('oxygen at 200m, \mu M','fontsize',14);
print -dpdf -painters fig1.pdf;

% difference plot
clear tmp
tmp(1:180,:)=o2obs(181:360,:);
tmp(181:360,:)=o2obs(1:180,:);
tmp(361,:)=tmp(1,:);
x=.5:360.5;
[xo,yo]=meshgrid(x,y);
o2obs2=griddata(xo,yo,tmp',xx,yy)';

load br.mat;
figure(2);
subplot1(2,2,'Gap',[.02 .02]);
for m=1:4
    del = o2mod(:,:,m) - o2obs2;
    subplot1(m);
    m_proj('robinson','clon',-150);
    m_pcolor(xx,yy,del');
    hold on;
    m_pcolor(xx-360,yy,del');
    hold off;
    shading flat;
    cmp=colormap(br);
    %m_grid('Xaxis','middle');
    m_grid('xticklabels','','yticklabels','');
    m_coast;
    caxis([-80 80]);
    
    a=o2mod(:,:,m);
    b=o2obs2;
    medo2(m)=nanmedian(a(:));
    corro2(m)=nancorr(a(:),b(:));
end
medo2(m+1)=nanmedian(o2obs2(:));

colorbartype([.1 .05 .8 .02],-80:4:80,40,[-80 80],cmp,0);
set(gca,'xtick',1:5:41);
set(gca,'xticklabel',{'-80' '-60' '-40' '-20' '0' ... 
    '+20' '+40' '+60' '+80'});
title('model - obs: oxygen, \mu M','fontsize',14);
print -dpdf -painters fig2.pdf;

% T comparison
x=-179.5:180.5;
[xo,yo]=meshgrid(x,y);
figure(3);
subplot1(3,2,'Gap',[0.02 0.02]);
subplot1(6);
m_proj('robinson','clon',-150);
m_pcolor(xo,yo,Tobs');
hold on;
m_pcolor(xo-360,yo,Tobs');
hold off;
shading flat;
caxis([0 25]);
cmp=colormap('jet');
m_grid('xticklabels','','yticklabels','');
m_coast;

for m=1:4
    subplot1(m);
    m_proj('robinson','clon',-150);
    m_pcolor(xx,yy,Tmod(:,:,m)');
    hold on;
    m_pcolor(xx-360,yy,Tmod(:,:,m)');
    hold off;
    shading flat;
    colormap('jet');
    m_grid('xticklabels','','yticklabels','');
    m_coast;
    caxis([0 25]);
end
subplot1(5);
set(gca,'visible','off');
colorbartype([.1 .2 .38 .02],0:1:25,26,[0 25],cmp,0);
set(gca,'xtick',1:5:26);
set(gca,'xticklabel',{'0' '5' '10' '15' '20' '25'});
title('T at 200m, deg C','fontsize',14);
print -dpdf -painters fig3.pdf;

% difference plot
clear tmp
tmp(1:180,:)=Tobs(181:360,:);
tmp(181:360,:)=Tobs(1:180,:);
tmp(361,:)=tmp(1,:);
x=.5:360.5;
[xo,yo]=meshgrid(x,y);
Tobs2=griddata(xo,yo,tmp',xx,yy)';

load br.mat;
figure(4);
subplot1(2,2,'Gap',[.02 .02]);
for m=1:4
    del = Tmod(:,:,m) - Tobs2;
    subplot1(m);
    m_proj('robinson','clon',-150);
    m_pcolor(xx,yy,del');
    hold on;
    m_pcolor(xx-360,yy,del');
    hold off;
    shading flat;
    cmp=colormap(br);
    %m_grid('Xaxis','middle');
    m_grid('xticklabels','','yticklabels','');
    m_coast;
    caxis([-4 4]);
    
    a=Tmod(:,:,m);
    b=Tobs2;
    medT(m)=nanmedian(a(:));
    corrT(m)=nancorr(a(:),b(:));
end
medo2(m+1)=nanmedian(o2obs2(:));

colorbartype([.1 .05 .8 .02],[-80:4:80]/20,40,[-4 4],cmp,0);
set(gca,'xtick',1:5:41);
set(gca,'xticklabel',{'-4' '-3' '-2' '-1' '0' ... 
    '+1' '+2' '+3' '+4'});
title('model - obs: T, deg C','fontsize',14);
print -dpdf -painters fig4.pdf;

% calculate MI
kh=O2sol(Sobs,Tobs)/0.21;
pO2obs=o2obs./kh;
khmod=O2sol(Smod,Tmod)/0.21;
pO2mod=o2mod./khmod;

% metabolic index 
E=[.6 .87 .74 .85 1.06 -.1 .68 .49 .56 .36 .23];
kb=8.6173303e-5;
taobs=Tobs+273.15;
tamod=Tmod+273.15;
E0(1)=prctile(E,25);
E0(2)=median(E);
E0(3)=prctile(E,75);

% Average size (B=3000 g) Atlantic Cod 
A=3.1e-14;
B=3e3;
ACobs=A*B^(-0.21)*pO2obs.*exp(E(2)./(kb*taobs));
ACmod=A*B^(-0.21)*pO2mod.*exp(E(2)./(kb*tamod));

% log o2 supply factor
Pobs = log(pO2obs);
Pmod = log(pO2mod);

% log o2 demand factor
Dobs = E0(2)./(kb.*taobs);

