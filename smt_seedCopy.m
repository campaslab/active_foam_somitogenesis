%% Copy original seed to surrounding boxes
function sdCp=smt_seedCopy(sd,bs)

% set copy_box of 8 surrounding squares
nFa=size(sd,1);
sdCp=zeros(9*nFa,2);
cpBx=bs*[0,0;1,1;1,0;1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];
    
% copy seed points to each surrounding square
for cbc=1:size(cpBx,1)
    sdCp(1+(cbc-1)*nFa:cbc*nFa,:)=sd+repmat(cpBx(cbc,:),nFa,1);        
end

end