%% Averaging intermediate vertices
% If intermediate vertices are very close, redistribute them to make almost
% uniform edge length

function [emd_n,rfn]=smt_edgeMidPtAvg(emd,edpc)

% compute edge length for each segment
eln=zeros(size(emd,1)-1,1);
for ii=1:size(emd,1)-1
    eln(ii)=norm(emd(ii+1,:)-emd(ii,:));
end

% erg is parametrization of each vertices between 0 to 1
rfn=floor(sum(eln)/edpc);
eln=eln/sum(eln);
erg=zeros(size(eln,1)+1,1);
for ii=2:size(erg)
    erg(ii)=sum(eln(1:ii-1));
end

% interpolate points based on parametric experssion
cnp=(1/(rfn+1))*(1:rfn);
nmd=zeros(size(cnp,2),2);
for ii=1:size(emd,1)-1
    for jj=1:size(nmd,1)
        nmd(jj,:)=nmd(jj,:)+...
            (emd(ii,:)+(cnp(jj)-erg(ii))/(erg(ii+1)-erg(ii))*...
            (emd(ii+1,:)-emd(ii,:)))*(heaviside(erg(ii+1)-cnp(jj))-...
            heaviside(erg(ii)-cnp(jj)));
    end
end    
emd_n=[emd(1,:);nmd;emd(end,:)];

end