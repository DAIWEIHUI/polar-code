function y=pencode(u,FZlookup,crc_size,bitreversedindices,F_kron_n) 


x = FZlookup;
switch crc_size
        case 0
            crc_code = [];
		case 4
			L=length(u);
            crc_gen=[1 0 0 1 1] ;        % CRC generator sequence
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);        % ����24λ��������2^24
            for i=1:L                               % ģ2��
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                   % CRC��
        case 6
			L=length(u);
            crc_gen=[1 1 0 0 0 0 1] ;        % CRC generator sequence
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);        % ����24λ��������2^24
            for i=1:L                               % ģ2��
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                   % CRC��
			
        case 8
            L=length(u);
            crc_gen=[1 0 0 0 0 0 1 1 1] ;        % CRC generator sequence
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);        % ����24λ��������2^24
            for i=1:L                               % ģ2��
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                   % CRC��
        case 11
            L=length(u);
            crc_gen=[1 1 1 0 0 0 1 0 0 0 0 1] ;        % CRC generator sequence
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);        % ����24λ��������2^24
            for i=1:L                               % ģ2��
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                   % CRC��
			
		case 12
		L=length(u);
		crc_gen=[1 1 0 0 0 0 0 0 0 1 1 1 1] ;        % CRC generator sequence
		left_shift=[1 zeros(1,crc_size)];
		a=conv(u,left_shift);        % ����24λ��������2^24
		for i=1:L                               % ģ2��
			if a(i)==1
				a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
			end
		end
		crc_code=a(L+1:L+crc_size);                   % CRC��


        case 16
            L=length(u);   
            crc_gen=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1];        % CRC generator sequence
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);               % ����16λ��������2^16
            for i=1:L                               % ģ2��
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                   % CRC��
            
        case 24
            L=length(u);   
            crc_gen=[1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1];        % CRC generator sequence
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);               % ����16λ��������2^16
            for i=1:L                               % ģ2��
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);                   % CRC��

        case 32
            L=length(u);
            crc_gen=[1 0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 0 1 1 0 1 1 1];        % CRC generator sequence
            left_shift=[1 zeros(1,crc_size)];
            a=conv(u,left_shift);        % ����16λ��������2^16
            for i=1:L                               % ģ2��
                if a(i)==1
                    a(i:i+crc_size)=xor(a(i:i+crc_size),crc_gen);
                end
            end
            crc_code=a(L+1:L+crc_size);   % CRC��
end
    
u = [u crc_code];

x (x == -1) = u; % -1's will get replaced by message bits below
x = x(bitreversedindices+1);
y = mod(x*F_kron_n,2);

end