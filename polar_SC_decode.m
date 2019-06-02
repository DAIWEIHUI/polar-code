function u_llr = polar_SC_decode(initial_llr,FZlookup,N,n)


    L = zeros(N,n+1);
    B = zeros(N,n+1);
    L(:,n+1) = initial_llr';
    
    for phi = 0:N-1
        [L,B] = updateL(n,phi,n,L,B);
        if FZlookup(phi+1) == 0
            B(phi+1,1) =  0;
        else
            if L(phi+1,1)<0
                B(phi+1,1) = 1;
            else
                B(phi+1,1) = 0;
            end
        end
        if mod(phi,2)==1
            [L,B] = updateB(n,phi,n,L,B);
        end
    end
    
    u_llr = L(:,1);
    
end


function [L,B] = updateL(lambda,phi,n,L,B)

    if lambda == 0
        return;
    end
    psi = floor(phi/2);
    if mod(phi,2)==0
        [L,B] = updateL(lambda-1,psi,n,L,B);
    end
    for omega=0:2^(n-lambda)-1
        if mod(phi,2)==0
            %do sth
            L(phi+omega*2^lambda+1,n+1-lambda) = fFunction(L(psi+2*omega*2^(lambda-1)+1,n+2-lambda),L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda));
        else
            %do sth
            if B(phi-1+omega*2^(lambda)+1,n+1-lambda) == 0
                L(phi+omega*2^(lambda)+1,n+1-lambda) = L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+L(psi+(2*omega)*2^(lambda-1)+1,n+2-lambda);
            else
                L(phi+omega*2^(lambda)+1,n+1-lambda) = L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)-L(psi+(2*omega)*2^(lambda-1)+1,n+2-lambda);
            end
        end
    end
end

function [L,B] = updateB(lambda,phi,n,L,B)
    
    psi = floor(phi/2);
    if mod(phi,2)~=0
        for omega = 0:2^(n-lambda)-1
            B(psi+2*omega*2^(lambda-1)+1,n+2-lambda) = xor(B(phi-1+omega*2^(lambda)+1,n+1-lambda),B(phi+omega*2^(lambda)+1,n+1-lambda));
            B(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda) = B(phi+omega*2^(lambda)+1,n+1-lambda);
        end
        if mod(psi,2)~=0
            [L,B] = updateB(lambda-1,psi,n,L,B);
        end
    end
end

function c = fFunction(a,b)
    c = sign(a).*sign(b).*min(abs(a),abs(b));
end



