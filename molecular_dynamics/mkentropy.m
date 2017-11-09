% compute quasiharmonic entropy
% supporting material for :
% Ovchinnikov, Louveau, Barton, Karplus and Chakraborty, submitted to eLife, 2017
%
run colours;

close all;

fname='entropy-ca';
fext='.mat' ;

ntraj=5 ;

tshc5=zeros(3,3,ntraj);
tslc5=zeros(3,3,ntraj);

for itraj=1:ntraj
% read entropy data domputed from 5 MD trajectories :
load([fname,num2str(itraj),fext]);

tshc=reshape(tsclass(:,1),3,3);
tshc5(:,:,itraj)=tshc ;

tslc=reshape(tsclass(:,2),3,3);
tslc5(:,:,itraj)=tslc ;


end

labels={'3BNC60/HC' 'CH103/HC' 'PGT121/HC'} ;
leg={'Mature','Intermediate','Germline'};

% average entropy over the five trajectories: 
tshc=mean(tshc5,3);
tshcstd=std(tshc5,0,3);

b=bar(tshc(end:-1:1,:)'+tshcstd(end:-1:1,:)', 0., 'linewidth',2); hold on % add 1 std
b=bar(tshc(end:-1:1,:)');

colormap( [ orange ; brown ; black ] );

ylabel('-TS^{cg}_{quasi}')
set(gca,'xticklabel',labels)

ylim([0 350]);
text(0.1,350,'A','fontsize',15);
%legend(leg,-1);

% generate figure
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc','quasi_hc_err.eps');

% same for light chain : 

labels={'3BNC60/LC' 'CH103/LC' 'PGT121/LC'} ;

figure;
%average entropy over the five trajectories
tslc=mean(tslc5,3);
tslcstd=std(tslc5,0,3);

b=bar(tslc(end:-1:1,:)'+tslcstd(end:-1:1,:)', 0., 'linewidth',2); hold on % add 1 std
b=bar(tslc(end:-1:1,:)');

colormap( [ orange ; brown ; black ] );

ylabel('-TS^{cg}_{quasi}')
set(gca,'xticklabel',labels)

ylim([0 350]);
text(0.1,350,'B','fontsize',15);

legend(leg(end:-1:1),2);legend boxoff;

% generate figure
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc','quasi_lc_err.eps');

