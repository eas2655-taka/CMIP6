%% validation
close all
clear all

% select reference level and time
z0=200;
t0=[1980 2000];
t1=[2080 2100];

% coordinates
x=-179.5:179.5;
y=-89.5:89.5;
xm=1:2:359;
ym=-89:2:89;
load validation/z33.mat;
z33=z;
time=[1850:2014]+0.5;
It=find(time>=t0(1)&time<=t0(2));
time=[2015:2099]+0.5;
It1=find(time>=t1(1)&time<=t1(2));

% open files 
dirs={'/Users/takamitsu3/Dropbox (GaTech)/Research/MI-CMIP6/WOA2018/' ... 
    '/Users/takamitsu3/Dropbox (GaTech)/Research/MI-CMIP6/CMIP6/'};
n0=[1 1 2 1];
model={'CanESM5' 'CNRM-ESM2-1' 'GFDL-CM4' 'IPSL-CM6A-LR'};
run={'r1i1p1f1' 'r1i1p1f2' 'r1i1p1f1' 'r1i1p1f1'};
run1={'r1i1p1f1' 'r4i1p1f2' 'r1i1p1f1' 'r1i1p1f1'};
sce={'ssp126' 'ssp245' 'ssp585'};
k=find(z33==z0);

for m=1:4
    % load model
    tmp=load([dirs{2},'/',model{m},'/o2_',model{m},'_historical_',run{m},'_2x2L33_ann_dft.mat']); 
    tmp2=mean(tmp.o2(:,:,:,It),4);
    for i=2:179
        for j=1:90
            if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
            end
        end
    end
    o2mod(:,:,m)=tmp2(:,:,k)*1e3;
    o2mod3d(:,:,:,m)=tmp2*1e3;
        
    clear tmp;
    tmp=load([dirs{2},'/',model{m},'/thetao_',model{m},'_historical_',run{m},'_2x2L33_ann_dft.mat']); 
    tmp2=mean(tmp.thetao(:,:,:,It),4);
    for i=2:179
        for j=1:90
            if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
            end
        end
    end
    Tmod(:,:,m)=tmp2(:,:,k);
    Tmod3d(:,:,:,m)=tmp2;

    clear tmp;
    tmp=load([dirs{2},'/',model{m},'/so_',model{m},'_historical_',run{m},'_2x2L33_ann_dft.mat']); 
    tmp2=mean(tmp.so(:,:,:,It),4);
    for i=2:179
        for j=1:90
            if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
            end
        end
    end
    Smod3d(:,:,:,m)=tmp2;
    Smod(:,:,m)=tmp2(:,:,k);
    
    for n=n0(m):3
        clear tmp;
        tmp=load([dirs{2},'/',model{m},'/thetao_',model{m},'_',sce{n},'_',run1{m},'_2x2L33_ann_dft.mat']); 
        tmp2=mean(tmp.thetao(:,:,:,It1),4);
        for i=2:179
            for j=1:90
                if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                    tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
                end
            end
        end
        Tmod1(:,:,n,m)=tmp2(:,:,k);
        Tmod3d1(:,:,:,n,m)=tmp2;
        
        clear tmp;
        tmp=load([dirs{2},'/',model{m},'/o2_',model{m},'_',sce{n},'_',run1{m},'_2x2L33_ann_dft.mat']); 
        tmp2=mean(tmp.o2(:,:,:,It1),4);
        for i=2:179
            for j=1:90
                if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                    tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
                end
            end
        end
        o2mod1(:,:,n,m)=tmp2(:,:,k)*1e3;
        o2mod3d1(:,:,:,n,m)=tmp2*1e3;

        clear tmp;
        tmp=load([dirs{2},'/',model{m},'/so_',model{m},'_',sce{n},'_',run1{m},'_2x2L33_ann_dft.mat']); 
        tmp2=mean(tmp.so(:,:,:,It1),4);
        for i=2:179
            for j=1:90
                if isnan(tmp2(i,j,1)) & ~isnan(tmp2(i+1,j,1)+tmp2(i-1,j,1))
                    tmp2(i,j,:)=.5*(tmp2(i+1,j,:)+tmp2(i-1,j,:));
                end
            end
        end
        Smod3d1(:,:,:,n,m)=tmp2;
        Smod1(:,:,n,m)=tmp2(:,:,k);
    end
end

% calculate MI
khmod=O2sol(Smod,Tmod)/0.21;
pO2mod=o2mod./khmod;
khmod1=O2sol(Smod1,Tmod1)/0.21;
pO2mod1=o2mod1./khmod1;

E=[.6 .87 .74 .85 1.06 -.1 .68 .49 .56 .36 .23];
kb=8.6173303e-5;
tamod=Tmod+273.15;
tamod1=Tmod1+273.15;
E0(1)=prctile(E,25);
E0(2)=median(E);
E0(3)=prctile(E,75);

% del (po2 factor)
do2 = permute(o2mod1,[1 2 4 3])-repmat(o2mod,[1 1 1 3]);
pO2mod1=permute(pO2mod1,[1 2 4 3]);
delo2 = real(log(pO2mod1)) - real(log(repmat(pO2mod,[1 1 1 3])));

% del (T factor)
dT = permute(Tmod1,[1 2 4 3])-repmat(Tmod,[1 1 1 3]);
tamod1=permute(tamod1,[1 2 4 3]);
tamod=repmat(tamod,[1 1 1 3]);
Tave = .5*(tamod + tamod1);
delT = -1./(kb*Tave.^2).*(tamod1 - tamod);

xm(end+1)=361;
[xx,yy]=meshgrid(xm,ym);
delo2(end+1,:,:,:)=delo2(1,:,:,:);
delT(end+1,:,:,:)=delT(1,:,:,:);

addpath ~/MATLAB/m_map
% O2 comparison, all ssps
figure(1);
subplot1(4,3,'Gap',[0.02 0.02]);
load br.mat;
cnt=0;
for m=1:4
    for n=1:3
        if m==3 & n== 1
            cnt=cnt+1;
            subplot1(cnt);
            set(gca,'visible','off');
        else
            cnt=cnt+1;
            subplot1(cnt);
            m_proj('robinson','clon',-150);
            m_pcolor(xx,yy,delo2(:,:,m,n)');
            hold on;
            m_pcolor(xx-360,yy,delo2(:,:,m,n)');
            hold off;
            shading flat;
            cmp=colormap(br);
            m_grid('xticklabels','','yticklabels','');
            m_coast;
            caxis([-.5 .5]);
        end
    end
end
colorbartype([.1 .05 .8 .02],-.5:.025:.5,41,[-.5 .5],cmp,0);
set(gca,'xtick',1:10:41);
set(gca,'xticklabel',{'-50%' '-25%' '0' '+25%' '+50%'});
title('fractional change in MI(pO2)','fontsize',14);
print -dpdf -painters fig5.pdf;

% T comparison, all ssps, E=0.3 (weak sensitivity)
figure(2);
subplot1(4,3,'Gap',[0.02 0.02]);
cnt=0;
for m=1:4
    for n=1:3
        if m==3 & n== 1
            cnt=cnt+1;
            subplot1(cnt);
            set(gca,'visible','off');
        else
            cnt=cnt+1;
            subplot1(cnt);
            m_proj('robinson','clon',-150);
            m_pcolor(xx,yy,.3*delT(:,:,m,n)');
            hold on;
            m_pcolor(xx-360,yy,.3*delT(:,:,m,n)');
            hold off;
            shading flat;
            cmp=colormap(br);
            m_grid('xticklabels','','yticklabels','');
            m_coast;
            caxis([-.5 .5]);
        end
    end
end
colorbartype([.1 .05 .8 .02],-.5:.025:.5,41,[-.5 .5],cmp,0);
set(gca,'xtick',1:10:41);
set(gca,'xticklabel',{'-50%' '-25%' '0' '+25%' '+50%'});
title('fractional change in MI(T, E=0.3)','fontsize',14);
print -dpdf -painters fig6.pdf;

% T comparison, all ssps, E=1.0 (strong sensitivity)
figure(3);
subplot1(4,3,'Gap',[0.02 0.02]);
cnt=0;
for m=1:4
    for n=1:3
        if m==3 & n== 1
            cnt=cnt+1;
            subplot1(cnt);
            set(gca,'visible','off');
        else
            cnt=cnt+1;
            subplot1(cnt);
            m_proj('robinson','clon',-150);
            m_pcolor(xx,yy,1*delT(:,:,m,n)');
            hold on;
            m_pcolor(xx-360,yy,1*delT(:,:,m,n)');
            hold off;
            shading flat;
            cmp=colormap(br);
            m_grid('xticklabels','','yticklabels','');
            m_coast;
            caxis([-.5 .5]);
        end
    end
end
colorbartype([.1 .05 .8 .02],-.5:.025:.5,41,[-.5 .5],cmp,0);
set(gca,'xtick',1:10:41);
set(gca,'xticklabel',{'-50%' '-25%' '0' '+25%' '+50%'});
title('fractional change in MI(T, E=1.0)','fontsize',14);
print -dpdf -painters fig7.pdf;

% T comparison, all ssps, E=0.6 (median sensitivity)
figure(4);
subplot1(4,3,'Gap',[0.02 0.02]);
cnt=0;
for m=1:4
    for n=1:3
        if m==3 & n== 1
            cnt=cnt+1;
            subplot1(cnt);
            set(gca,'visible','off');
        else
            cnt=cnt+1;
            subplot1(cnt);
            m_proj('robinson','clon',-150);
            m_pcolor(xx,yy,.6*delT(:,:,m,n)');
            hold on;
            m_pcolor(xx-360,yy,.6*delT(:,:,m,n)');
            hold off;
            shading flat;
            cmp=colormap(br);
            m_grid('xticklabels','','yticklabels','');
            m_coast;
            caxis([-.5 .5]);
        end
    end
end
colorbartype([.1 .05 .8 .02],-.5:.025:.5,41,[-.5 .5],cmp,0);
set(gca,'xtick',1:10:41);
set(gca,'xticklabel',{'-50%' '-25%' '0' '+25%' '+50%'});
title('fractional change in MI(T, E=0.6)','fontsize',14);
print -dpdf -painters fig8.pdf;


% net MI comparison, all ssps, E=0.6 (median sensitivity)
figure(5);
subplot1(4,3,'Gap',[0.02 0.02]);
cnt=0;
net = .6*delT + delo2;
for m=1:4
    for n=1:3
        if m==3 & n== 1
            cnt=cnt+1;
            subplot1(cnt);
            set(gca,'visible','off');
        else
            cnt=cnt+1;
            subplot1(cnt);
            m_proj('robinson','clon',-150);
            m_pcolor(xx,yy,net(:,:,m,n)');
            hold on;
            m_pcolor(xx-360,yy,net(:,:,m,n)');
            hold off;
            shading flat;
            cmp=colormap(br);
            m_grid('xticklabels','','yticklabels','');
            m_coast;
            caxis([-.5 .5]);
        end
    end
end
colorbartype([.1 .05 .8 .02],-.5:.025:.5,41,[-.5 .5],cmp,0);
set(gca,'xtick',1:10:41);
set(gca,'xticklabel',{'-50%' '-25%' '0' '+25%' '+50%'});
title('fractional change in MI(net, E=0.6)','fontsize',14);
print -dpdf -painters fig9.pdf;

% zonal mean summary plots
figure(6);
subplot1(3,3,'Gap',[.02 .02]);
delo2x = squeeze(nanmean(delo2,1))*1e2;
delTx = squeeze(nanmean(delT,1))*1e2;
netx=delo2x+.6*delTx;
delTx(:,3,:)=NaN;
for m=1:3
    subplot1(3*m-2);
    plot(ym,delo2x(:,:,m)); %legend(p,model);
    axis([-70 70 -100 60]);
    grid on;

    subplot1(3*m-1);
    plot(ym,delTx(:,:,m));
    axis([-70 70 -100 60]);
    grid on;

    subplot1(3*m);
    p=plot(ym,netx(:,:,m));
    axis([-70 70 -100 60]);
    grid on;
end
legend(p,model);
orient landscape
print -dpdf -painters fig10.pdf;

% o2-t ratio
figure(7);
R = do2./dT;
R(abs(dT)<.25)=NaN;
Rx=squeeze(nanmean(R,1));

