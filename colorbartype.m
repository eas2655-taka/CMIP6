function [h0]=colorbartype(pos,lev,lev_pas,clip_lev,hpal,orien);

% 30/10/97- NS: option "orien" pour positionner la barre verticalement
%               ou horizontalement--> 1 ou 0

h0 = axes('Position',pos,'Box','on');
ind_lev=find( lev>=(clip_lev(1)-eps) & (lev<=clip_lev(2)+eps) );
tab_lev=[lev(ind_lev);lev(ind_lev)];

% horizontale par defaut
if (nargin < 6), orien = 0; end 

if orien == 0,

    pcolor(tab_lev);
    set(h0,'xticklabel', '','yticklabel','');
    abs_bar=0:1/(length(ind_lev)-1):1;
    for j=1:lev_pas:length(ind_lev)
        s=num2str(lev(ind_lev(j)));
%        text(abs_bar(j),-.5,s,'Fontsize', 12, 'HorizontalAlignment', ...
%                              'Center', 'Units', 'normalized');
    end

else

    pcolor(tab_lev');
    set(h0,'xticklabel', '','yticklabel','');
    ord_bar=0:1/(length(ind_lev)-1):1;
    for j=1:lev_pas:length(ind_lev)
        s=num2str(lev(ind_lev(j)));
%        text(1.5,ord_bar(j),s,'Fontsize', 12,'Units', 'normalized');
    end

end
colormap(hpal);
caxis(clip_lev);

