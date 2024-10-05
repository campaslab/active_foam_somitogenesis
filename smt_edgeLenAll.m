function eLnc=smt_edgeLenAll(edge,gm_p)    

eLnc=zeros(gm_p.nEd,1);
for edc=1:gm_p.nEd
    eLnc(edc)=sum(edge{4}{edc});
end

end