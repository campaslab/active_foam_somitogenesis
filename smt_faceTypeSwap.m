function face_n=smt_faceTypeSwap(edge,face,ftp1,ftp2,gm_p,rg)

face_n=face;

% Swap one face from type 1 to type 2.
faId=find(face_n{1}(:,1)==ftp1);
fneId=cell(size(faId,1),1);
[fneCnt,fneCnt2]=deal(zeros(size(faId,1),1)); 

for fac=1:size(fneId,1)
    fneId{fac}=unique(edge{1}(abs(face{3}{faId(fac)}),3:4));
    fneId{fac}=fneId{fac}(fneId{fac}~=faId(fac));
    fneCnt2(fac)=sum(fneId{fac}==0);
    fneId{fac}=fneId{fac}(fneId{fac}~=0);
    fneId{fac}=face_n{1}(fneId{fac},1);
    fneCnt(fac)=sum(fneId{fac}==ftp2);
end

faSwapId=find(fneCnt>1 & fneCnt2==0);
faSwapId=faSwapId(randperm(size(faSwapId,1)));
faSwapId=faId(faSwapId(1));    
face_n{1}(faSwapId,1)=ftp2;

% Swap one face from type 2 to type 1.
faId=find(face_n{1}(:,1)==ftp2);
fneId=cell(size(faId,1),1);
[fneCnt,fneCnt2]=deal(zeros(size(faId,1),1)); 

for fac=1:size(fneId,1)
    fneId{fac}=unique(edge{1}(abs(face_n{3}{faId(fac)}),3:4));
    fneId{fac}=fneId{fac}(fneId{fac}~=faId(fac));
    fneCnt2(fac)=sum(fneId{fac}==0);
    fneId{fac}=fneId{fac}(fneId{fac}~=0);
    fneId{fac}=face_n{1}(fneId{fac},1);
    fneCnt(fac)=sum(fneId{fac}==ftp1);
end

faSwapId=find(fneCnt>1 & fneCnt2==0);
faSwapId=faSwapId(randperm(size(faSwapId,1)));
faSwapId=faId(faSwapId(1));

face_n{1}(faSwapId,1)=ftp1;

% delete isolate cells
face_n=smt_isolatedFaceRsgn(edge,face_n,gm_p,rg);
% faId=find(face_n{1}(:,1)==ftp1);
% for fac=1:size(faId,1)
%     fneId=unique(edge{1}(abs(face_n{3}{faId(fac)}),3:4));
%     fneId=fneId(fneId~=0 & fneId~=faId(fac));
%     fneId=face_n{1}(fneId,1);
%     if sum(fneId==ftp1)==0
%         face_n{1}(faId(fac))=ftp2;
%     end
% end
% 
% faId=find(face_n{1}(:,1)==ftp2);
% for fac=1:size(faId,1)
%     fneId=unique(edge{1}(abs(face_n{3}{faId(fac)}),3:4));
%     fneId=fneId(fneId~=0 & fneId~=faId(fac));
%     fneId=face_n{1}(fneId,1);
%     if sum(fneId==ftp2)==0
%         face_n{1}(faId(fac))=ftp1;
%     end
% end  

end