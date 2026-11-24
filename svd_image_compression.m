%% ========================================================================
%  IMAGE COMPRESSION USING SINGULAR VALUE DECOMPOSITION (SVD)
%  ------------------------------------------------------------------------
%  A grayscale image is a matrix A (m x n). SVD factorises it as
%
%        A = U * S * V'
%
%  where the diagonal of S holds the singular values in decreasing order.
%  Keeping only the k largest singular values gives the best possible
%  rank-k approximation of the image (Eckart-Young theorem):
%
%        A_k = U(:,1:k) * S(1:k,1:k) * V(:,1:k)'
%
%  A_k needs only k*(m + n + 1) numbers instead of m*n, so the image is
%  compressed while keeping most of its visual quality.
%
%  How to run: open this file in MATLAB and press Run (or type
%  svd_image_compression in the Command Window). All figures and results
%  are also saved in the "output" folder.
%  Works in MATLAB (no toolbox needed) and in GNU Octave.
%  ========================================================================

clc; clear; close all;

%% ---------------- Settings (change these if you like) -------------------
imageFile   = 'sample.png';          % any JPG/PNG image in this folder
kValues     = [5 10 20 50 100];      % numbers of singular values to keep
targetPSNR  = 30;                    % dB - commonly used "good quality" level
outputDir   = 'output';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% ---------------- Step 1: Load the image and convert to grayscale -------
img = imread(imageFile);

if size(img, 3) == 3                 % colour image -> grayscale
    if exist('rgb2gray')             %#ok<EXIST> % available in MATLAB
        grayImg = rgb2gray(img);
    else                             % same formula rgb2gray uses
        grayImg = uint8(0.2989 * double(img(:,:,1)) + ...
                        0.5870 * double(img(:,:,2)) + ...
                        0.1140 * double(img(:,:,3)));
    end
else
    grayImg = img;                   % already grayscale
end

fig1 = figure('Name', 'Step 1 - Original and Grayscale Image', 'Position', [100 100 1000 420], 'PaperPositionMode', 'auto');
subplot(1, 2, 1); imshow(img);     title('Original Image');
subplot(1, 2, 2); imshow(grayImg); title('Grayscale Image');
print(fig1, fullfile(outputDir, 'fig1_original_grayscale.png'), '-dpng', '-r150');

%% ---------------- Step 2: Represent the image as a matrix ---------------
A = double(grayImg);                 % pixel values 0-255 as a double matrix
[m, n] = size(A);
r = rank(A);

fprintf('=====================================================\n');
fprintf(' IMAGE COMPRESSION USING SVD\n');
fprintf('=====================================================\n');
fprintf('Image file            : %s\n', imageFile);
fprintf('Matrix size (m x n)   : %d x %d\n', m, n);
fprintf('Total pixel values    : %d\n', m * n);
fprintf('Rank of image matrix  : %d\n', r);
fprintf('Sample of matrix A (top-left 5 x 5 pixel values):\n');
disp(A(1:5, 1:5));

%% ---------------- Step 3: Apply SVD -------------------------------------
[U, S, V] = svd(A);
sigma = diag(S);                     % singular values, largest first

fprintf('Size of U: %d x %d,  S: %d x %d,  V: %d x %d\n', ...
        size(U,1), size(U,2), size(S,1), size(S,2), size(V,1), size(V,2));
fprintf('Number of singular values : %d\n', numel(sigma));
fprintf('Largest 10 singular values:\n');
fprintf('  %10.2f\n', sigma(1:10));

% Link with eigenvalues: singular values of A are the square roots of the
% eigenvalues of A'*A.
lambda = sort(eig(A' * A), 'descend');
lambda = lambda(1:numel(sigma));
eigCheck = max(abs(sqrt(max(lambda, 0)) - sigma) ./ sigma(1));
fprintf('Check sigma = sqrt(eig(A''*A)): max relative difference = %.2e\n', eigCheck);

% Check that U*S*V' rebuilds A exactly
fprintf('Check U*S*V'' = A: max reconstruction difference = %.2e\n\n', ...
        max(max(abs(U * S * V' - A))));

% Energy (information) captured by the first k singular values
energy    = sigma .^ 2;
cumEnergy = cumsum(energy) / sum(energy) * 100;

fig2 = figure('Name', 'Step 3 - Singular Values', 'Position', [100 100 1100 450], 'PaperPositionMode', 'auto');
subplot(1, 2, 1);
semilogy(sigma, 'b-', 'LineWidth', 1.5); grid on;
xlabel('Index i'); ylabel('Singular value \sigma_i (log scale)');
title('Singular Values of the Image');
subplot(1, 2, 2);
plot(cumEnergy, 'r-', 'LineWidth', 1.5); grid on;
xlabel('Number of singular values k'); ylabel('Energy retained (%)');
title('Cumulative Energy'); ylim([min(cumEnergy) - 1, 100.5]);
print(fig2, fullfile(outputDir, 'fig2_singular_values.png'), '-dpng', '-r150');

%% ---------------- Step 4: Reconstruct with different k -----------------
nK        = numel(kValues);
relError  = zeros(1, nK);   % relative Frobenius error (%)
mseVals   = zeros(1, nK);   % mean squared error
psnrVals  = zeros(1, nK);   % peak signal-to-noise ratio (dB)
storage   = zeros(1, nK);   % numbers stored for the compressed image
ratio     = zeros(1, nK);   % compression ratio
saved     = zeros(1, nK);   % space saved (%)
energyK   = zeros(1, nK);   % energy retained (%)
rankK     = zeros(1, nK);   % rank of the compressed matrix
compressed = cell(1, nK);

originalStorage = m * n;

for i = 1:nK
    k  = kValues(i);
    Ak = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';    % rank-k approximation
    Ak = min(max(Ak, 0), 255);                     % keep valid pixel range
    compressed{i} = uint8(Ak);
    imwrite(compressed{i}, fullfile(outputDir, sprintf('compressed_k%03d.png', k)));

    % ------------ Step 5: Error and storage for this k -------------------
    diffMat     = A - Ak;
    relError(i) = norm(diffMat, 'fro') / norm(A, 'fro') * 100;
    mseVals(i)  = mean(diffMat(:) .^ 2);
    psnrVals(i) = 10 * log10(255^2 / mseVals(i));
    storage(i)  = k * (m + n + 1);                 % U_k + V_k + k singular values
    ratio(i)    = originalStorage / storage(i);
    saved(i)    = (1 - storage(i) / originalStorage) * 100;
    energyK(i)  = cumEnergy(k);
    rankK(i)    = rank(U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)');
end

%% ---------------- Step 5: Error analysis and storage comparison ---------
fprintf('-----------------------------------------------------------------------------------------\n');
fprintf('%5s %6s %12s %10s %10s %10s %10s %10s %10s\n', 'k', 'Rank', 'Storage', ...
        'Ratio', 'Saved %', 'Energy %', 'Error %', 'MSE', 'PSNR dB');
fprintf('-----------------------------------------------------------------------------------------\n');
fprintf('%5s %6d %12d %10s %10s %10s %10s %10s %10s\n', 'Full', r, originalStorage, ...
        '1.00', '0.00', '100.00', '0.00', '0.00', 'Inf');
for i = 1:nK
    fprintf('%5d %6d %12d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n', kValues(i), ...
            rankK(i), storage(i), ratio(i), saved(i), energyK(i), relError(i), ...
            mseVals(i), psnrVals(i));
end
fprintf('-----------------------------------------------------------------------------------------\n');

% Error curve for every k, using the Eckart-Young result:
% ||A - A_k||_F^2 = sum of the squares of the discarded singular values
allK        = 1:numel(sigma);
tailEnergy  = sum(energy) - cumsum(energy);
tailEnergy(tailEnergy < 0) = 0;
relErrAll   = sqrt(tailEnergy) / norm(A, 'fro') * 100;
psnrAll     = 10 * log10(255^2 ./ (tailEnergy / (m * n)));
ratioAll    = originalStorage ./ (allK * (m + n + 1));

% Smallest k that reaches the target quality
kGood = find(psnrAll >= targetPSNR, 1);
breakEven = floor(originalStorage / (m + n + 1));
fprintf('Smallest k with PSNR >= %d dB : k = %d  (PSNR %.2f dB, compression ratio %.2f, energy %.2f%%)\n', ...
        targetPSNR, kGood, psnrAll(kGood), ratioAll(kGood), cumEnergy(kGood));
fprintf('Compression saves space only while k < %d (break-even point)\n\n', breakEven);

fig4 = figure('Name', 'Step 5 - Reconstruction Error', 'Position', [100 100 1100 450], 'PaperPositionMode', 'auto');
plotK = 1:min(200, numel(sigma));
subplot(1, 2, 1);
plot(plotK, relErrAll(plotK), 'b-', 'LineWidth', 1.5); hold on;
plot(kValues, relError, 'ro', 'MarkerFaceColor', 'r');
grid on; xlabel('Number of singular values k'); ylabel('Relative error (%)');
title('Reconstruction Error vs k'); legend('All k', 'Chosen k');
subplot(1, 2, 2);
plot(plotK, psnrAll(plotK), 'b-', 'LineWidth', 1.5); hold on;
plot(kValues, psnrVals, 'ro', 'MarkerFaceColor', 'r');
plot([plotK(1) plotK(end)], [targetPSNR targetPSNR], 'k--');
grid on; xlabel('Number of singular values k'); ylabel('PSNR (dB)');
title('Image Quality (PSNR) vs k');
legend('All k', 'Chosen k', sprintf('%d dB target', targetPSNR), 'Location', 'southeast');
print(fig4, fullfile(outputDir, 'fig4_error_analysis.png'), '-dpng', '-r150');

fig5 = figure('Name', 'Step 5 - Storage Comparison', 'Position', [100 100 1100 450], 'PaperPositionMode', 'auto');
subplot(1, 2, 1);
bar([originalStorage, storage] / 1000);
set(gca, 'XTickLabel', [{'Orig'}, arrayfun(@(k) sprintf('%d', k), kValues, 'UniformOutput', false)]);
xlabel('k (singular values kept)'); ylabel('Values stored (thousands)'); title('Storage Required'); grid on;
subplot(1, 2, 2);
bar(ratio);
set(gca, 'XTickLabel', arrayfun(@(k) sprintf('%d', k), kValues, 'UniformOutput', false));
xlabel('k (singular values kept)'); ylabel('Compression ratio'); title('Compression Ratio'); grid on;
print(fig5, fullfile(outputDir, 'fig5_storage_comparison.png'), '-dpng', '-r150');

%% ---------------- Step 6: Display original and compressed images --------
fig3 = figure('Name', 'Step 6 - Original vs Compressed Images', 'Position', [50 50 1200 720], 'PaperPositionMode', 'auto');
subplot(2, 3, 1); imshow(grayImg);
title({'Original', sprintf('rank %d', r)});
for i = 1:min(nK, 5)
    subplot(2, 3, i + 1); imshow(compressed{i});
    title({sprintf('k = %d', kValues(i)), ...
           sprintf('Ratio %.1f : 1,  PSNR %.1f dB', ratio(i), psnrVals(i))});
end
print(fig3, fullfile(outputDir, 'fig3_compressed_images.png'), '-dpng', '-r150');

%% ---------------- Save results to a CSV file ----------------------------
fid = fopen(fullfile(outputDir, 'results.csv'), 'w');
fprintf(fid, 'k,Rank,ValuesStored,CompressionRatio,SpaceSavedPercent,EnergyRetainedPercent,RelativeErrorPercent,MSE,PSNR_dB\n');
for i = 1:nK
    fprintf(fid, '%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', kValues(i), rankK(i), ...
            storage(i), ratio(i), saved(i), energyK(i), relError(i), mseVals(i), psnrVals(i));
end
fclose(fid);

fprintf('All figures, compressed images and results.csv saved in the "%s" folder.\n', outputDir);
