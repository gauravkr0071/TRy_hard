




%% Neighborhood Matrix Generation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dated: 08/02/2017

% Copyright (C) Raghav Dev, All Right Reserved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vector,index, Window]=R_ij_M(image,i,j,M)
x=(i-M):(i+M);
% fprintf('x is %d \n',size(x));
y=(j-M):(j+M);
% fprintf('y is %d \n',size(y));
neighborhood_length = (2*M+1)^2;
% fprintf('neighborhood length is %d \n',neighborhood_length);
Window=image(x,y);
% fprintf('size of window is %d \n',size(Window));
vector=reshape(Window,[1,neighborhood_length]);
index=combvec(x,y);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%






