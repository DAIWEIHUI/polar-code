function [u_llr, c_llr] = polar_SCAN_decode(y_llr,iter_num,FZlookup,N)
    %初始化PCparams.L 和 PCparams.B

    
    n = log2(N);


    plus_infinity = 1000;
    L = zeros(N,n+1);%left message
    B = zeros(N,n+1);%right message
    L(:,n+1) = y_llr';%initial L
    B(FZlookup==0,1) = plus_infinity;%initial B
    
    
    %主循环
    for ii = 1:iter_num
        for phi = 0:N-1
            [L,B] = updateLLRMap(n,phi,n,L,B);
            if mod(phi,2)~=0
                [L,B] = updateBitMap(n,phi,n,L,B);   
            end
        end
    end

    %输出最终的左信息u_llr和右信息c_llr
    u_llr = L(:,1)'+B(:,1)';
    
    c_llr = B(:,n+1)';
    
    
end


function [L,B] = updateLLRMap(lambda,phi,n,L,B)

    if lambda == 0
        return;
    end
    psi = floor(phi/2);
    if mod(phi,2)==0
        [L,B] = updateLLRMap(lambda-1,psi,n,L,B);
    end
    for omega=0:2^(n-lambda)-1
        if mod(phi,2)==0
            %do sth
            L(phi+omega*2^lambda+1,n+1-lambda) = fFunction(L(psi+2*omega*2^(lambda-1)+1,n+2-lambda),L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+B(phi+1+omega*2^lambda+1,n+1-lambda));
        else
            %do sth
            L(phi+omega*2^lambda+1,n+1-lambda) = L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+fFunction(L(psi+2*omega*2^(lambda-1)+1,n+2-lambda),B(phi-1+omega*2^lambda+1,n+1-lambda));
        end
    end
end


function [L,B] = updateBitMap(lambda,phi,n,L,B)

    
    psi = floor(phi/2);
    if mod(phi,2)~=0
        for omega = 0:2^(n-lambda)-1
            B(psi+2*omega*2^(lambda-1)+1,n+2-lambda) = fFunction(B(phi-1+omega*2^(lambda)+1,n+1-lambda),B(phi+omega*2^(lambda)+1,n+1-lambda)+L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda));
            B(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda) = B(phi+omega*2^(lambda)+1,n+1-lambda) + fFunction(B(phi-1+omega*2^(lambda)+1,n+1-lambda),L(psi+2*omega*2^(lambda-1)+1,n+2-lambda));
        end
        if mod(psi,2)~=0
            [L,B] = updateBitMap(lambda-1,psi,n,L,B);
        end
    end
end

function c = fFunction(a,b)
    c = sign(a)*sign(b)*min(abs(a),abs(b));
end