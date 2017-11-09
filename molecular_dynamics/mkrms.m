% plot rmsfs
% supporting material for :
% Ovchinnikov, Louveau, Barton, Karplus and Chakraborty, submitted to eLife, 2017

close all;

allnames=[ { '3bnc60t' '3bnc60at', '3bnc60glt'}
           {'ch103t', 'ch103-i3.2t', 'ch103ucat'}
           {'pgt121t', '3h109lt', 'gl121t'} ];

for names=allnames'

% names
% continue

% legend names
 legnames={'Mature','Intermediate','Germline'};

 flags={'-cg-hc', '-cg-lc' };
 lflgs={'heavy', 'light' };
 run colours;
 clrs={ [black] [brown] [orange]};

 aligns={'malign.mat' 'malign-lc.mat'}; % heavy and light chain sequence alignments, respectively
 xlims=[145,120]; % axis limits
 ylims=[6,6];

 lw=1.5;
%
 for jj=1:length(flags)
  leg={};
  hleg=[];
  flag=char(flags(jj)); % file name flag
  lflg=char(lflgs(jj));
%
  align=char(aligns(jj));
  load(align); % load seq. aligmnents
  armsf=zeros(size(maind));
% indicate domains
  figure('position',[100 300 900 300]); hold on; box on;
  showdom; % draw domain limits
% plot RMSFs
  for ii=1:3
    armsf(:)=NaN; % do not forget to reinitialize
    name=char(names(ii));
    clr=cell2mat(clrs(ii)) ; % color array

% load rmsf data for this antibody
    load(['./',name,'-rmsf',flag,'.mat']);
    n=length(flucall);
    resnum=[1:n];

% low pass filter to make plots clearer :
    d=3; % use a frame of 3 adjacent residues
    fs =smooth2(resnum,flucall,d);
    fsp=smooth2(resnum,flucall+flucstd,d);
    fsm=smooth2(resnum,flucall-flucstd,d);
%
% put all rmsf into one variable ( armsf ) for plotting, using the sequence alignment read above
% find the correct antibody sequence index in maind
    ind=find(ismember(files,name(1:end-1))); % end-1 to drop 't' et the end
    for k=1:length(armsf)
     i3=maind(k,ind);
     if (i3>0 & i3<=n)
      armsf(k,1)=fs(i3);
      armsf(k,2)=fsp(i3);
      armsf(k,3)=fsm(i3);
     end
    end

  h=plot(armsf(:,1), 'color', clr, 'linewidth',lw) ;hold on ; box on ;
    plot(armsf(:,2),'--', 'color', clr,'linewidth',1) ;hold on ; box on ;
    plot(armsf(:,3),'--', 'color', clr,'linewidth',1) ;hold on ; box on ;
%
    leg = [ leg {[char(legnames(ii))]} ];
    hleg=[hleg h]; %legend handles

% mark residues that are mutated in comparison with the germline/common ancestor
% using sequence alignment read above
    ms=6;
    if ( strcmp( name(1:6),'3bnc60') )
     seq0=ma(3,:); % 3bnc gl
     panels='AD';
    elseif ( strcmp( name(1:5),'ch103') )
     seq0=ma(5,:); % for ch103 uca
     panels='BE';
    else
     seq0=ma(9,:) ; % pgt
     panels='CF';
    end

    if (ii==1 | ii==2)
     seqm=ma(ind,:);
     mind = find ( (seq0~=seqm).*(seq0~='-').*(seqm~='-') ) ;
     for k=mind
      c=cell2mat(clrs(ii));
      plot(k,armsf(k,1), 'o', 'color',c,'markersize',ms,'markerfacecolor',c)
% omit labels since they are hard to see
%   text(k,armsf(k,ii), [seq0(k),seqm(k)])
     end
    end
%===
  end
% generate figure
  xlabel('\it residue', 'fontsize',14);
  ylabel('$\it RMSF(\AA)$', 'interpreter','latex', 'fontsize',14);
% draw subfigure letter
  text(-8, 6, ['\it',panels(jj),')'], 'fontsize',14.5);
  xlim([0 xlims(jj)]);
  ylim([0 ylims(jj)]);
  legend(hleg,leg,-1,'fontsize',14);
  set(gcf,'paperpositionmode','auto');
  print(gcf, '-depsc2', [char(names(1)),'rmsf',flag,'.eps']);
 end

end