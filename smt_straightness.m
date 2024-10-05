function strt=smt_straightness(vrtx,edge,face,gm_p,ftp1,ftp2)

efTp=zeros(gm_p.nEd,2);
efTp(edge{1}(:,3)~=0,1)=face{1}(edge{1}(edge{1}(:,3)~=0,3),1);
efTp(edge{1}(:,4)~=0,2)=face{1}(edge{1}(edge{1}(:,4)~=0,4),1);
efTp=sort(efTp,2);

bndLen=0;
ftp=sort([ftp1,ftp2]);
bndId=find(efTp(:,1)==ftp(1) & efTp(:,2)==ftp(2));
for edc=1:size(bndId,1)
    bndLen=bndLen+sum(edge{4}{bndId(edc)});
end

vfTp=zeros(gm_p.nVr,3);
for ii=1:3
    vfTp(vrtx(:,ii+3)~=0,ii)=face{1}(vrtx(vrtx(:,ii+3)~=0,ii+3),1);
end
vfTp=sort(vfTp,2);

bndVrId=find(vfTp(:,1)==0 & vfTp(:,2)==ftp(1) & vfTp(:,3)==ftp(2));
if size(bndVrId,1)~=2
    error('There must be only two end vertices!');
end

vrCrd=vrtx(bndVrId,7:8);
strt=norm(vrCrd(2,:)-vrCrd(1,:))/bndLen;

end