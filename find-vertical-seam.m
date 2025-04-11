% Step 1: Read the Image
I = imread('cat.jpg'); 
I_gray = rgb2gray(I); 

% Step 3: Add Noise to the Image
noise_density = 0.1;
I_noisy = imnoise(I_gray, 'salt & pepper', noise_density);
I_gray = I_noisy;

% Step 2: Compute Edges
dx = [1, 0, -1];
dy = [1; 0; -1];

Gx = imfilter(double(I_gray), dx, 'replicate');
Gy = imfilter(double(I_gray), dy, 'replicate');


% Step 3: Calculate Energy using L1 Norm
energy_map = abs(Gx) + abs(Gy);

% Step 4: Normalize Energy Map for Display
energy_map = mat2gray(energy_map);

% Step 5: Initialize Cumulative Energy Table
[rows, cols] = size(energy_map);
cumulative_energy = zeros(rows, cols);
cumulative_energy(1, :) = energy_map(1, :);


% Step 4: Compute Cumulative Energy
for r = 2:rows
    for c = 1:cols
        % Handle boundary conditions
        left = max(c - 1, 1);
        right = min(c + 1, cols);
        
        % Find the minimum energy path to the current pixel
        cumulative_energy(r, c) = energy_map(r, c) + ...
            min(cumulative_energy(r-1, left:right));
    end
end

% Step 5: Backtrack to Find Optimal Seam
% Start from the last row, find the minimum energy
[~, min_col] = min(cumulative_energy(end, :));
seam = zeros(rows, 1);
seam(rows) = min_col;

for r = rows-1:-1:1
    c = seam(r+1);
    left = max(c - 1, 1);
    right = min(c + 1, cols);
    [~, idx] = min(cumulative_energy(r, left:right));
    seam(r) = left + idx - 1; % Adjust index
end


% Step 8: Visualize the Seam
output_image = cat(3, energy_map, energy_map, energy_map);
for r = 1:rows
    % Mark the seam in red
    output_image(r, seam(r), 1) = 255; % Red channel
    output_image(r, seam(r), 2) = 255; % Green channel
    output_image(r, seam(r), 3) = 0;   % Blue channel
end

figure;
imagesc(energy_map);
colormap(jet);
colorbar;
title('Energy Map');
axis image;

% Step 9: Display Results
figure;
subplot(1, 2, 1);
imshow(I_noisy);
title('Noisy Image');

subplot(1, 2, 2);
imshow(output_image);
title('Optimal Vertical Seam');
