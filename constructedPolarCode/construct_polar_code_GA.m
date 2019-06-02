%使用GA方法构造polar code
%polar encoder and SC, BP , SCAN decoder for awgn channel
%panzhipeng
function construct_polar_code_GA(N,sigma)
    n = ceil(log2(N)); 
    NN = 2^n;
    if(NN~=N)
        fprintf('The num N must be the power of 2!');
        return;
    end
    
    file_name = sprintf('PolarCode_block_length_%d_sigma_%.2f_method_GA.txt',N,sigma);
    fid = fopen(file_name,'w+');
    bitreversedindices = zeros(1,N);
    for index = 1 : N
        bitreversedindices(index) = bin2dec(wrev(dec2bin(index-1,n)));
    end
    initialize_phi();% construct table of phi(x) and inv(-log(phi(x)));
    

    mean_llr = 2/sigma^2;
    channels = mean_llr*ones(1, N);
    
    for i=1:n
        c1 = channels(1:2:N);
        c2 = channels(2:2:N);
        channels = [phi_x_inv(1 - (1 - phi_x_table(c1)).*(1 - phi_x_table(c2)))', c1 + c2];
    end
    channels = channels(bitreversedindices+1);

    [~,indices] = sort(channels,'descend');
    for ii = 1:length(channels)
        fprintf(fid,'%d\r\n',indices(ii));
    end
    fclose(fid);
end


function [ x ] = phi_x_inv( y )

global minus_log_phi_inv_table;
global min_minus_log_phi;
global max_minus_log_phi;
global increment_minus_log_phi;
minus_log_phi = -log(y);
minus_log_phi = max(minus_log_phi, min_minus_log_phi);
minus_log_phi = min(minus_log_phi, max_minus_log_phi);
minus_log_phi_index = round((minus_log_phi - min_minus_log_phi)/increment_minus_log_phi - 0.499) + 1;
x = minus_log_phi_inv_table(minus_log_phi_index);

end


function [ y ] = phi_x_table( x )

global phi_x_table;
global min_x;
global max_x;
global increment_x;
x = max(min_x, x);
x = min(max_x, x);
x_index = round((x - min_x)/increment_x) + 1;
y = phi_x_table(x_index);

end

function  initialize_phi(x_increment)

persistent phi_initialized;
global minus_log_phi_inv_table;
global min_minus_log_phi;
global max_minus_log_phi;
global increment_minus_log_phi;
global phi_x_table;
global min_x;
global max_x;
global increment_x;

if isempty(phi_initialized)
    
    phi_initialized = 1;
    if nargin < 1
        x_increment = 10^(-4);
    end
    x_vec = 0: x_increment: 400;
    
    
    phi_vec = (x_vec <  10) .* exp(-0.4527 * x_vec.^(0.86) + 0.0218);
    phi_vec = (x_vec >= 10) .* sqrt(pi./(x_vec+0.0001)) .* (1 - 1.4286./(x_vec + 0.0001)) .* exp(-x_vec/4) + phi_vec;
    phi_vec = min(phi_vec, 1);
    min_minus_log_phi = 0;
    max_minus_log_phi = 100;
    increment_minus_log_phi = 10^(-3);
    minus_log_phi_inv_table = zeros(ceil((max_minus_log_phi - min_minus_log_phi)/increment_minus_log_phi) + 1, 1);
    
    min_x = 0;
    max_x = 100;
    increment_x = 0.01;
    index = 1;
    for x = min_x: increment_x: max_x + increment_x
        if (x < 10)
            phi_x_table(index) = exp(-0.4527 * x.^(0.86) + 0.0218);
        else
            phi_x_table(index) =sqrt(pi./(x)) .* (1 - 1.4286./(x)) .* exp(-x/4);
        end
        index = index + 1;
    end
    
    
    for x_index = 1 : length(x_vec)
        x_value  = x_vec(x_index);
        phi_value = phi_vec(x_index);
        minus_log_phi = -log(phi_value);
        if minus_log_phi < max_minus_log_phi + increment_minus_log_phi
            minus_log_phi_index = ceil((minus_log_phi - min_minus_log_phi)/increment_minus_log_phi) + 1;
            minus_log_phi_inv_table(minus_log_phi_index) = x_value;
        end
    end
    
end

end

