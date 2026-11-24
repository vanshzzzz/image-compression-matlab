# Image Compression Using Singular Value Decomposition (SVD)

This MATLAB mini project compresses a grayscale image by keeping only the k largest
singular values of its matrix. It shows the original and compressed images, the
reconstruction error and the storage savings.

## Files

| File | What it is |
|---|---|
| `svd_image_compression.m` | The complete MATLAB program (all 6 steps) |
| `sample.png` | Sample photo (coffee cup, 400 × 600, CC0 / no copyright restrictions) |
| `output/` | Results from a test run: 5 figures, 5 compressed images, `results.csv`, `console_output.txt` |

## How to run

### MATLAB (desktop)
1. Extract the zip.
2. Open MATLAB. In the **Current Folder** panel, go to the extracted `svd_image_compression` folder.
3. Double-click `svd_image_compression.m` to open it.
4. Click **Run** (green ▶ button in the Editor tab), or type this in the Command Window:
   ```
   svd_image_compression
   ```
5. Five figure windows open and the results table prints in the Command Window.
   Everything is also saved in the `output` folder.

### MATLAB Online (no installation)
1. Go to **matlab.mathworks.com** and sign in with your college email.
2. Click **Upload** and upload `svd_image_compression.m` and `sample.png`.
3. Open the `.m` file and click **Run**.

### GNU Octave (free alternative)
Open Octave, `cd` to the folder, then run `svd_image_compression`.

## Using your own image
Copy any `.jpg` or `.png` into the folder and change the first setting in the script:
```
imageFile = 'myphoto.jpg';
```
You can also change which k values are tested:
```
kValues = [5 10 20 50 100];
```

## What the program does

| Step | Action |
|---|---|
| 1 | Loads the image with `imread()` and converts it to grayscale |
| 2 | Converts it to a matrix `A` and prints its size and rank |
| 3 | Applies `[U,S,V] = svd(A)`, plots the singular values and cumulative energy, and checks that σ = √eig(AᵀA) |
| 4 | Rebuilds the image with k = 5, 10, 20, 50, 100 singular values |
| 5 | Calculates relative error, MSE, PSNR, storage, compression ratio and space saved |
| 6 | Displays the original and compressed images with `imshow()` |

## Results on the sample image (400 × 600, rank 400)

| k | Values stored | Compression ratio | Space saved | Relative error | PSNR |
|---|---|---|---|---|---|
| Original | 240,000 | 1.00 | 0 % | 0 % | ∞ |
| 5 | 5,005 | 47.95 | 97.91 % | 22.10 % | 19.74 dB |
| 10 | 10,010 | 23.98 | 95.83 % | 17.00 % | 22.02 dB |
| 20 | 20,020 | 11.99 | 91.66 % | 13.24 % | 24.20 dB |
| 50 | 50,050 | 4.80 | 79.15 % | 9.40 % | 27.17 dB |
| 100 | 100,100 | 2.40 | 58.29 % | 5.99 % | 31.09 dB |

- **k = 86** is the smallest k with PSNR ≥ 30 dB (good quality), at about 2.8× compression.
- Compression saves space only while **k < 239**: m·n / (m + n + 1) = 240000 / 1001.
- Storage is counted in **numbers stored**. Each k keeps one column of U (m values),
  one column of V (n values) and one singular value, so k(m + n + 1) in total.

## Key formulas

```
A   = U * S * V'                               (SVD)
A_k = U(:,1:k) * S(1:k,1:k) * V(:,1:k)'        (rank-k approximation)
Storage(k)        = k * (m + n + 1)
Compression ratio = (m * n) / Storage(k)
Relative error    = ||A - A_k||_F / ||A||_F * 100
MSE  = mean((A - A_k).^2)
PSNR = 10 * log10(255^2 / MSE)
Energy retained   = sum(sigma(1:k).^2) / sum(sigma.^2) * 100
```
