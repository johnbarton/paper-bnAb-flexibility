% draw domain limitss

ind=[mafwr ; macdr ] ; % assume that we have read malign.mat ; framework and cdr region limits

legends={ 'FWR1', 'FWR2', 'FWR3', 'FWR4', 'CDR1', 'CDR2', 'CDR3' };

run colours
colors=repmat(white,8,1);

str=ind(:,1); %starting indices
stp=ind(:,2); %stopping indices

% note that the above domain definitions are with respect to 3bnc sequence
% need to convert to unwrapped/aligned sequence :
seqind=3 ;% 3bnc60 gl
inda=zeros(size(ind)) ;
for i=1:length(ind)
 for j=1:2
  inda(i,j)=find(maind(:,seqind)==ind(i,j));
 end
end
%
str=inda(:,1); %starting indices
stp=inda(:,2); %stopping indices

%plot the lines

l=length(str);
yl=6;

cols=colormap;
cinc=length(cols)/l;

gray=[0.8 0.8 0.8];
white=[1 1 1];
red=[1 0 0];

for i=1:l
 x=linspace(str(i), stp(i),2);
 h=rectangle('position', [ x(1) 0 x(2)-x(1) yl ], 'facecolor', colors(i,:), 'edgecolor',[0 0 0]);
% write text labels
 t=text( 0.5 * ( x(2) + x(1) ) + 0 , 4.75 , char(legends(i)), 'fontsize',12,'rotation', 90,'fontweight','bold');
%
end

