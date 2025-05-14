function Main
Para=1;
Repeats=32; % Number of 3D structures to generate
%%%

% Load the input SEM image
A=imread(['Raw/' 'DN.tif']);%%%%%
if size(A, 3) == 3
   A = rgb2gray(A);  % Convert to grayscale if RGB
end
[rows, cols] = size(A);
min_dimension_value = min(rows, cols);
pl = input("Enter the pixel length in microns: ");
if min_dimension_value >= 300
    updated_pixel_length = (300/min_dimension_value) * pl;
else
    updated_pixel_length = (min_dimension_value/300) * pl;
end
disp("Updated Pixel Length: " + updated_pixel_length + " microns");
user_response = input("Enter 'y' to continue: ", 's');
if strcmpi(user_response, 'y')
    disp("Continuing with the code...");
else
    disp("Code execution stopped.");
    return;
end
%%%

% Perform reconstruction in parallel or serial

if Para==1
parfor II=1:Repeats
    Calc(II);
end
end
if Para==0
for II=1:Repeats
    Calc(II);
end
end
end

function Calc(II)
close all
opt.show=0;
opt.Steps=200;
rng(randi(1e8));
A=imread(['Raw/' 'DN.tif']);%%%%%%%%%
if numel(size(A))==3; A=rgb2gray(A); end
A=cropBinaryImage(A);

A=imresize(A,[430,430]);  %%%% Resize the initial SEM image (adjusting pixel size accordingly to the desired piel size) %%%%

A = imgaussfilt(A, 1);  % Smooth image to reduce noise
B = imbinarize(A, 'adaptive', 'sensitivity', 0.61);  % Adaptive binarization
B = medfilt2(B, [3, 3]);  % Median filter to clean binary image

S=[size(A,1),size(A,2),3]; % 3D size of output volume
opt.F0=0;
opt.S=S;
opt.B0=B;
X = find(B);  % Run optimization to find best parameters
C = gen(X, S);  % Generate 3D structure using parameters

% Save central slice side-by-side with input binary
G=cat(2,C(:,:,2),B);
imwrite(uint8(G.*255),['Rec/' num2str(II) '.tif']);
save(['X/' num2str(II) '.mat'],'X')
close all;

    function X=find(A)
        % Bayesian optimization of parameters X1 to X9
        opt.iter=0;
        S=[size(A,1),size(A,2),5];
        F0=feature(A);
        opt.F0=F0;
        opt.S=S;
        opt.MinE=1e6;
        x1 = optimizableVariable('x1',[0,1],'Type','real');
        x2 = optimizableVariable('x2',[0,1],'Type','real');
        x3 = optimizableVariable('x3',[0,1],'Type','real');
        x4 = optimizableVariable('x4',[0,1],'Type','real');
        x5 = optimizableVariable('x5',[0,1],'Type','real');
        x6 = optimizableVariable('x6',[0,1],'Type','real');
        x7 = optimizableVariable('x7',[0,1],'Type','real');
        x8 = optimizableVariable('x8',[0,1],'Type','real');
        x9 = optimizableVariable('x9',[0,1],'Type','real');
        MAXEVAL=opt.Steps;
        results = bayesopt(@Err,[x1,x2,x3,x4,x5,x6,x7,x8,x9],...
            'ExplorationRatio',.99,...
            'AcquisitionFunctionName','expected-improvement',...
            'MaxObjectiveEvaluations',MAXEVAL,...
            'MaxTime',150,...
            'GPActiveSetSize',50,'PlotFcn',...
            {@plotMinObjective,@plotAcquisitionFunction});
        X=table2array(results.XAtMinObjective);
    end

    function E=Err(X)
        % Objective function: mean absolute error between feature vectors
        opt.iter=opt.iter+1;
        if istable(X); X=table2array(X); end
        AT=gen(X,opt.S);
        F=feature(AT);
        try
            t1=opt.F0+.01;
            t2=F+.01;
            E=mean(abs(t1-t2));  % Mean absolute error (MAE) as optimization objective
            if opt.show
                subplot(2,3,1); imagesc(opt.B0); axis tight equal;
                subplot(2,3,2); imagesc(AT(:,:,ceil(size(AT,3)/2))); axis tight equal;
            end
            if opt.MinE>E
                opt.MinE=E;
                opt.BestF=F;
                opt.BestX=X;
                if opt.show
                    subplot(2,3,3);
                    imagesc(AT(:,:,ceil(size(AT,3)/2))); axis tight equal;
                end
            end
            if opt.show
                subplot(2,3,4); cla; plot(opt.F0); hold on; title(num2str(E)); plot(opt.BestF,'color','g'); title(num2str(E)); drawnow;
                subplot(2,3,5); hold on; scatter(opt.iter,opt.MinE,'k');
                subplot(2,3,6); cla; bar([opt.BestX ; X]'); ylim([0,1]);
            end
        catch
            E=1e6; % Penalize failed evaluations
        end
    end
end

function A = gen(X, S) % Generate 3D porous structure from parameter vector X
if nargin < 2
    S = [256, 256, 6];
end
Margin = 10;
Margin2 = 10;
MaxSigma = max(S) / 25;
X1 = X(1);
X2 = X(2);
X3 = X(3);
X4 = X(4);
X5 = X(5);
X6 = X(6);
X7 = X(7);
X8 = X(8);
X9 = X(9);

% AM1: Gaussian-filtered noise field
AM1 = rand([S + 2 * Margin]);
AM1 = imgaussfilt3(AM1, [X1 * MaxSigma + 1.5, X1 * MaxSigma + 1.5, X1 * MaxSigma + 1.5],'padding','symmetric');
AM1 = AM1(Margin + 1:S(1) + Margin, Margin + 1:S(2) + Margin, Margin + 1:S(3) + Margin);

% AM2: Distance map from rare events
AM2 = rand([S + 2 * Margin2]);
AM2 = bwdist((AM2 < 10^(-(X2 * 3 + 2.5))));
AM2 = AM2(Margin2 + 1:S(1) + Margin2, Margin2 + 1:S(2) + Margin2, Margin2 + 1:S(3) + Margin2);
A1 = AM1 < quantile(AM1(:), X3 * 0.9 + 0.05);

% Construct combined maps B1, B2, B3
B1 = bwdist(1 - A1) - bwdist(A1);
B1 = normalize(B1);
B2 = AM2;
B2 = normalize(max(B2(:))-B2);
AM3 = rand([S + 2 * Margin]);
AM3 = bwdist((AM3 < 10^(-(X9 * 3 + 2.5))));
AM3 = AM3(Margin + 1:S(1) + Margin, Margin + 1:S(2) + Margin, Margin + 1:S(3) + Margin);
B3 = AM3;
B3 = normalize(B3);

% Combine and threshold the blended volume
B = B1*X4+ B2*X8+B3*X7;
B = imgaussfilt3(B, [X5 * MaxSigma, X5 * MaxSigma, X5 * MaxSigma],'padding','symmetric') - B;
B = imgaussfilt3(B, [.5,.5,.5],'padding','symmetric');
A = B < quantile(B(:), 1 - (X6 * 0.7 + 0.05));
A=medfilt3(A,[3,3,3]);
end

function [F]=feature(A)
% Extracts multi-scale histogram features for similarity comparison
% Compute histograms of distances for both phases at 3 zoom levels
if numel(size(A))==3; A=A(:,:,ceil(size(A,3)/2)); end
B1=bwdist(A); Q1=B1(:); Q1(Q1==0)=[];
B2=bwdist(~A); Q2=B2(:); Q2(Q2==0)=[];
F1=[fliplr(histcounts(Q1,1:.2:25)) , histcounts(Q2,1:.2:25)]./numel(A);
A=imzoom(A,1.5);
B1=bwdist(A); Q1=B1(:); Q1(Q1==0)=[];
B2=bwdist(~A); Q2=B2(:); Q2(Q2==0)=[];
F2=[fliplr(histcounts(Q1,1:.2:25)) , histcounts(Q2,1:.2:25)]./numel(A);
A=imzoom(A,1.5);
B1=bwdist(A); Q1=B1(:); Q1(Q1==0)=[];
B2=bwdist(~A); Q2=B2(:); Q2(Q2==0)=[];
F3=[fliplr(histcounts(Q1,1:.2:50)) , histcounts(Q2,1:.2:50)]./numel(A);
F=[F1 F2 F3];
end

function B=imzoom(A,Rat)
% Zoom image by padding and resizing
m=ceil((Rat-1)./2.*size(A));
B=imresize(A,size(A)+2.*m,'nearest');
B=B(m(1)+1:end-m(1),m(2)+1:end-m(2));
end
function croppedImage = cropBinaryImage(image)
[rows, cols] = size(image);
if rows > cols
croppedImage = image(1:cols, :);
elseif cols > rows
croppedImage = image(:, 1:rows);
else
croppedImage = image;
end
end