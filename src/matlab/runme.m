%   Copyright (C) 2015, Gabriele Facciolo

% test 1: compare SGM (MGM=1), basic MGM (MGM=2), and MGM with 4 messages (MGM=4)
figure(1)
L=imread('data/imL.png');
R=imread('data/imR.png');
dmax=16;
P1=10;
P2=20;
DIR=8;

i=1;
for MGM=[1,2,4] 
   [dmap,t] = stereomatch_MGM(L, R, dmax, DIR, P1, P2, MGM);
   subplot(1,3,i); imagesc(dmap); axis image;
   title([' MGM:' num2str(MGM) ' t:' num2str(t)])
   i=i+1;
end


% test 2: compare SGM, MGM, and MGM with truncated linear potential
figure(2)
L=imread('data/fountain23-imL.png');
R=imread('data/fountain23-imR.png');
dmax=143;
P1=4;
P2=16;
DIR=4; 
MGM=2;
i=1;

% baseline SGM 
[dmap,t] = stereomatch_MGM(L, R, dmax, DIR, P1, P2, 1, 0);
subplot(1,3,i); imagesc(dmap); axis image;
title(['TLP:' num2str(0) ' MGM:' num2str(MGM) ' t:' num2str(t)])
i=i+1;
% TEST MGM and MGM with truncated linear potential
for TLP=[0,1] 
      if(TLP==1) 
         PP2=P1*15;
      else
         PP2=P1*4;
      end
      [dmap,t] = stereomatch_MGM(L, R, dmax, DIR, P1, PP2, MGM, TLP);
      subplot(1,3,i); imagesc(dmap); axis image;
      title(['TLP:' num2str(TLP) ' MGM:' num2str(MGM) ' t:' num2str(t)])
      i=i+1;
end
