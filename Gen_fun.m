

% Define the reconstruction parameters (X1 to X9)
% These control noise filtering, thresholding, weighting, smoothing, etc.

X=[0.9187237835497712,0.3858329603184352,0.7922278198405153,0.15925835080401848,0.46428111000860467,0.23866660081994182,0.9714286072553279,0.8145359319775307,0.40271441892950705],



% Define the size of the generated 3D structure
%%%% These dimensions match the resized 2D SEM image in the main script%%%%
S=[430,430,430];


% Choose a random slice index to extract from the 3D volume
rn=randi([0 200]);

% Generate the 3D structure
A = gen(X, S);
ccc=A*255;
imwrite(ccc(:,:,rn), '2D.tif');

imwrite(ccc(:,:,1), '3D.tif');


for i=2:size(ccc,3)
    imwrite(ccc(:,:,i), '3D.tif', 'WriteMode', 'append');
end

% --- 3D Generation Function ---
function A = gen(X, S)
if nargin < 2
    S = [256, 256, 6];
end
Margin = 10;
Margin2 = 10;
MaxSigma = max(S) / 25; % Max smoothing sigma scaled to volume size

X1 = X(1);
X2 = X(2);
X3 = X(3);
X4 = X(4);
X5 = X(5);
X6 = X(6);
X7 = X(7);
X8 = X(8);
X9 = X(9);

% First Gaussian noise field (base structure generation)
AM1 = rand([S + 2 * Margin]);
AM1 = imgaussfilt3(AM1, [X1 * MaxSigma + 1.5, X1 * MaxSigma + 1.5, X1 * MaxSigma + 1.5],'padding','symmetric');
AM1 = AM1(Margin + 1:S(1) + Margin, Margin + 1:S(2) + Margin, Margin + 1:S(3) + Margin);

% Second field: rare-event binary mask distance map
AM2 = rand([S + 2 * Margin2]);
AM2 = bwdist((AM2 < 10^(-(X2 * 3 + 2.5))));
AM2 = AM2(Margin2 + 1:S(1) + Margin2, Margin2 + 1:S(2) + Margin2, Margin2 + 1:S(3) + Margin2);

% Binary thresholding of AM1 to form structure A1
A1 = AM1 < quantile(AM1(:), X3 * 0.9 + 0.05);

% Generate B1: distance map difference across void/solid
B1 = bwdist(1 - A1) - bwdist(A1);
B1 = normalize(B1);

% Generate B2: normalized inverse of AM2
B2 = AM2;
B2 = normalize(max(B2(:))-B2);

% Third noise map (similar to AM2, using different sparsity level)
AM3 = rand([S + 2 * Margin]);
AM3 = bwdist((AM3 < 10^(-(X9 * 3 + 2.5))));
AM3 = AM3(Margin + 1:S(1) + Margin, Margin + 1:S(2) + Margin, Margin + 1:S(3) + Margin);
B3 = AM3;
B3 = normalize(B3);

% Combine maps using learned weights
B = B1*X4+ B2*X8+B3*X7;

% Apply 3D Gaussian smoothing and normalize again
B = imgaussfilt3(B, [X5 * MaxSigma, X5 * MaxSigma, X5 * MaxSigma],'padding','symmetric') - B;
B = imgaussfilt3(B, [.5,.5,.5],'padding','symmetric');

% Global threshold to generate binary volume
A = B < quantile(B(:), 1 - (X6 * 0.7 + 0.05));

% Apply 3D median filtering to clean isolated pixels
A=medfilt3(A,[3,3,3]);
end

