%% Plot all cells

function smt_plot(edge,face,rg,gm_p)

% color scheme for each region
col=[{'r'},{'g'},{'b'},{'c'},{'y'},{'m'}];
colSz=size(col,2);

% draw individual cells with patch function
figure;
for fac=1:gm_p.nFa
    fvrCrd=smt_faceVrtx(face{3}{fac},edge{2});
    ftpInd=mod(face{1}(fac,rg.fi(1)),colSz)...
        +colSz*(mod(face{1}(fac,rg.fi(1)),colSz)==0);    
    patch(fvrCrd(:,1),fvrCrd(:,2),col{ftpInd});
    alpha(0.5);
    hold on;
end

% set axis 
axis([-2 gm_p.bs+2 -2 gm_p.wth+2]);
pbaspect([gm_p.bs+4 gm_p.wth+4 1]);

end