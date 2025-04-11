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


% Step 7: Compute Cumulative Energy
for c = 2:cols
    for r = 1:rows
        % Handle boundary conditions
        up = max(r - 1, 1);
        down = min(r + 1, rows);
        
        % Find the minimum energy path to the current pixel
        cumulative_energy(r, c) = energy_map(r, c) + ...
            min(cumulative_energy(up:down, c-1));
    end
end

% Step 8: Backtrack to Find Optimal Seam
% Start from the last column, find the minimum energy
[~, min_row] = min(cumulative_energy(:, end));
seam = zeros(1, cols);
seam(end) = min_row;

for c = cols-1:-1:1
    r = seam(c+1);
    up = max(r - 1, 1);
    down = min(r + 1, rows);
    [~, idx] = min(cumulative_energy(up:down, c));
    seam(c) = up + idx - 1; % Adjust index
end

% Step 9: Visualize the Seam
output_image = cat(3, energy_map, energy_map, energy_map);
for c = 1:cols
    % Mark the seam in red
    output_image(seam(c), c, 1) = 255; % Red channel
    output_image(seam(c), c, 2) = 255; % Green channel
    output_image(seam(c), c, 3) = 0;   % Blue channel
end

% Step 10: Display Results
figure;
imagesc(energy_map);
colormap(jet);
colorbar;
title('Energy Map');
axis image;

figure;
subplot(1, 2, 1);
imshow(I_noisy);
title('Noisy Image');

subplot(1, 2, 2);
imshow(output_image);
title('Optimal Horizontal Seam');
