
close all;
clear;
clc;

%% Salt and Pepper Noise removal Using type 2 fuzzy system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dated: 08/02/2017
% Copyright (C) Raghav Dev, Vikas Singh, Narendra K. Dhar;
% All Right Reserved.
% Last Modified on: 16/02/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is to filter GRAY image noised with salt 
% and pepper noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% time=tic;


input_image=imread('Lenna.png');
figure;
imshow(input_image);
title('Input image');

im_gray=rgb2gray(input_image);
% im_gray=input_image;
im_gray_1=im2double(im_gray); % Will be easy to calculate PSNR in the end

figure;
imshow(im_gray_1);
title('Input gray image');

%% Inialization Of Parameters
Noise_density = 0.50;
im_noised=imnoise(im_gray_1,'salt & pepper',Noise_density);
figure;
imshow(im_noised);
title(sprintf('Input noisy image with %d noise density',Noise_density));
[p,q]=size(im_noised); 

% Inialising Denoised image %% Image size is being increased by 8, to take
% care of edges.

im_denoised=0.63*ones(p+16,q+16); 
im_denoised(9:p+8,9:q+8)=im_noised; %

% M=1;  % M is half window size 
epsilon=0.0001;
count0=0;
count1=0;
count2=0;
count3=0;

% N_init=4;   % for low noise image N can be low
% N_init=8; % for higher noise image N should be high
M=1;  % M is half window size %% Inialization
N_init = 8;   % for low noise image N can be low

im_denoised_pixels = zeros(p+16,q+16);
im_noised_pixels = zeros(p+16,q+16);
%%
time_elapsed_per_itteration = zeros(1,10);
tic
e=1;
for z=1:2 % For loop for Iterations
    im_denoised_pixels = zeros(p+16,q+16);
    im_noised_pixels = zeros(p+16,q+16);
    %     c=dm;  % To terminate the loop (not necessary to use)
    for j=9:q+8 % To scan rows
        for i=9:p+8 % To scan col.
            M = 1;  % M is half window size %% Inialization
            N_init = 6;   % for low noise image N can be low
%             
            S_max = 2; % upper bound of 'M'
            N = N_init; % Inialisation of number of good pixel required to evaluate the pixel value of a currupted pixel
            
            while (im_denoised(i,j)==0)||(im_denoised(i,j)==1)
%                 M_Matrix(e)=M;
            %    e=e+1;
                [R_ij_M_matrix,index] = R_ij_M(im_denoised,i,j,M); % Vector containing pixel around possible currupted pixel
                
                [T_min,T_max,T_min_max,PI,H,sigma,average_mu] = membership_type_2(R_ij_M_matrix);
%                 [T_min,T_max,T_min_max,PI,H]=membership_type_2_1(R_ij_M_matrix);
                lenght_R_ij_M_matrix = length(R_ij_M_matrix);
                ave_PI = sum(PI)/H;
                % T_Threshold=(T_min+T_max+T_min_max)/3;
                T_Threshold = T_max;
                o = 1;
                %G = zeros(1,lenght_R_ij_M_matrix);
                %w = G;
                if ave_PI(H)>T_Threshold
                    count0=count0+1;
                    break
                elseif sigma==epsilon
                    % im_denoised(i,j) = average_mu(1);
                    im_denoised_pixels(i,j) = average_mu(1);
                    im_noised_pixels(i,j) = im_denoised(i,j);
                    count1=count1+1;
                    break
                end
                G=zeros(1,N);
                for x=1:lenght_R_ij_M_matrix
                     if  ave_PI(x)>=T_Threshold % (ave_PI(x)<T_Threshold)&&(ave_PI(x)>T_min)
                         count2=count2+1;
                         G(o)=R_ij_M_matrix(x);
                         o=o+1;
%                         w(x)=1/((((i-index(1,x))^2)+((j-index(2,x))^2))^2.5);
%                         % w(x)=1/(0.5*M*M);
                     elseif  (R_ij_M_matrix(x)~=0)&&(R_ij_M_matrix(x)~=1)
                    % if  (R_ij_M_matrix(x)~=0)&&(R_ij_M_matrix(x)~=1)
                        count3=count3+1;
                        G(o)=R_ij_M_matrix(x);
                        o=o+1;
                       % w(x)=1/((((i-index(1,x))^2)+((j-index(2,x))^2))^2.5);
                        % w(x)=1/(0.5*M*M);
                    end
                    if o==N+1
                        break
                    end
                end
                
                neta=length(find(G));
                
                if (neta<N) && (M< S_max)
                    M=M+1;                    
                    continue
                elseif (neta<N) && (M== S_max)
                    N=N-1;
                    if N<1
                        S_max=S_max+1;
                        N=1;
                    end
                    continue
                end
                for k=1:(length(G)/2)
                    mean(k)=K_middle_mean_w(k,G);
                end
                mean_G=((sum(mean))/length(G));
                var_G=2.5*abs(G-mean_G);
                var_G=max(var_G);
                if var_G<=.01
                    var_G=.01;
                end
                w=gaussmf(G,[var_G,mean_G]);
                W=sum(w);
                weighted_G=w*G';
                % im_denoised(i,j)=weighted_G/W;
                im_denoised_pixels(i,j) = weighted_G/W;
                im_noised_pixels(i,j) = im_denoised(i,j);
                break
            end
        end
    end
    time_elapsed_per_itteration(z)=toc;
    im_denoised=(im_denoised-im_noised_pixels)+im_denoised_pixels;
end
time_elapsed=toc;

%%

im_denoised=im_denoised(9:p+8,9:q+8);

im_denoised_1 = im2uint8(im_denoised);

im_denoised_1=int16(im_denoised_1);

im_gray=int16(im_gray);

PSNR=10*log10((255*255)/((1/((p-10)*(q-10)))*sum(sum((im_denoised_1(6:p-5,6:q-5)-im_gray(6:p-5,6:q-5)).^2))));

% PSNR=10*log10(1/((1/((p-14)*(q-14)))*sum(sum((im_gray(8:p-7,8:q-7)-im_denoised(8:p-7,8:q-7)).^2))));
% PSNR=10*log10(1/((1/((p)*(q)))*sum(sum((im_gray-im_denoised).^2))));

% PSNR_char=num2str(PSNR);
fprintf('PSNR is %d\n',PSNR);

figure;
imshow(im_denoised);
title(sprintf('Denoised image with %d noise density',Noise_density));
% title(PSNR_char);

