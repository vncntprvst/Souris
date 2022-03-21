fid=[whiskerData.fid];
wid=[whiskerData.wid];
labels=[whiskerData.label];
folY=[whiskerData.follicle_y];

fIdx=1:find(fid==251,1)-1;
wid=wid(fIdx);
labels=labels(fIdx);
folY=folY(fIdx);

figure; hold on 
plot(folY(wid==0))
plot(folY(wid==1))

figure; hold on 
plot(folY(labels==0))
plot(folY(labels==1))


