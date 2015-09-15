function [dmap,T] = stereomatch_MGM(imgleft, imgright, disparity, NDIR, P1, P2, MGM, VTYPE)

% stereomatch_MGM - Stereo with More Global Matching (MGM)
%
%   [dmap,T] = stereomatch_MGM(imgleft, imgright, disparity, [NDIR(8)], [P1(8)], [P2(32)], [MGM(2)])
%
%   Fills cost volume (using AD) and calls MGM to perform the stereo matching.
%   http://dev.ipol.im/~facciolo/mgm
%   
%   INPUTS:
%   imgleft/imgright: images of the stereo pair
%   disparity: maximum disparity (final range is [0, disparity])
%   NDIR : number of propagation directions as in SGM (default: 8, i.e. N,S,E,W)
%   P1/P2: parameters of interaction potential (default: 8/32)
%   MGM  : number of neighbours used for the propagation, 1 yield SGM (default: 2)
%   VTYPE:  1 yield SGM (default: 2)
%
%   OUTPUTS:
%   dmap: is the disparity map computed by MGM 
%   T   : processing time
%
%   Copyright (C) 2015, Gabriele Facciolo


	if nargin < 4
		NDIR = 8;
	end
	if nargin < 5
		P1 = 8.0;
	end
	if nargin < 6
		P2 = 32.0;
	end
	if nargin < 7
		MGM = 2;
    end
	if nargin < 8
		VTYPE = 0;
	end

   % Set Parameters
   D   = uint16(disparity)+1;              % number of disparities
   heightL = size( imgleft, 1 );
   widthL  = size( imgleft, 2 );
   
   % Initialization
   pcost = ones( heightL, widthL, D, 'single' )*255*255;
   
   % Calculate pixel cost
   for Dc = 1 : D
      maxL = widthL + 1 - Dc; 
      pcost(:, Dc : widthL, Dc ) = mean(imabsdiff( imgright( :, 1 : maxL), imgleft( :, Dc : widthL) ),3);
   end

   % adapt edge weights (TODO)
   w=ones(heightL,widthL,8);
   % w(:,2:widthL,1)   = mean( imabsdiff(imgleft(:,2:widthL,:), imgleft(:,1:widthL-1,:)),3);
   % w(:,1:widthL-1,2) = mean( imabsdiff(imgleft(:,1:widthL-1,:), imgleft(:,2:widthL,:)),3);
   % w(1:heightL-1,:,3)= mean( imabsdiff(imgleft(1:heightL-1,:,:), imgleft(2:heightL,:,:)),3);
   % w(2:heightL,:,4)  = mean( imabsdiff(imgleft(2:heightL,:,:), imgleft(1:heightL-1,:,:)),3);
   % w...,5)...
   % w = exp(-w.*w/64)*4;
   
   % Call the MGM_wrapper
   tic;
   dmap = MGM_wrapper(pcost, NDIR, P1, P2, MGM, VTYPE, w);
   T=toc;
   return
