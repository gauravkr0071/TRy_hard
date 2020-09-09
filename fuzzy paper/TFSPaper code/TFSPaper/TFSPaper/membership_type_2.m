%% Neighborhood Matrix Generation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dated: 08/02/2017

% Copyright (C) Raghav Dev, Vikas Singh, and Narendra K. Dhar
% All Right Reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 



function [T_min,T_max,T_min_max,PI,H,sigma,average_mu]=membership_type_2(vector)
p=numel(vector);
H=(p+1)/2;

% Because K=1 results mdian of set which cound be possibly 1 or 0 for high
% level of noise, hence it it better to start with K=2; hence lenght of
% mean vector will be one less that half lenght of vector.

mu=zeros(1,H);
for q=1:H;
    mu(1,q)=K_middle_mean(q,vector);
end


average_mu=((sum(mu))/H)*ones(1,p);
% minimum_mu=min(mu); % Need to justify our choice!!!!!!!

lembda= 1.5*abs(vector-average_mu).^2;

%
% lembda= abs(vector-average_mu);
% lembda= (vector-average_mu).^2; % uNCOMMENT FOR USING SQUERE OF DIFFERENCE
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PI is matrix which will save all the membership values of 
% each element of any vector. one element of vector will have
% H number of membership values corresponding to 'H' number of
% mean i. e. mu.

% PI=zeros(H,p);

mu_Mat=repmat(mu,p,1)';
lembda_vector=repmat(vector,H,1);

sigma=K_middle_mean(H,lembda);

if sigma < 0.0001
    sigma=0.0001;
end

PI= exp(-0.5*((lembda_vector-mu_Mat)/sigma).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% max_PI=max(PI);
% min_Pi=min(PI);
T_max=max(max(PI));

T_min_max=min(max(PI));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "max of min" can be used as threshold for membership value 
% to decide which is good and which is bad pixel

T_min=max(min(PI)); %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end










