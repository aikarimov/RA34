%Fstr = 'ddJFFF + dJJFF + 2dJFJF + JdJFF + JJJF';
Fstr = 'dddJFFFF + ddJJFFF + 2ddJFJFF + 3ddJFFJF + dJdJFFF + dJJJFF + 3dJJFJF + 3dJFdJFF + 3dJFJJF + JddJFFF + JdJJFF + 2JdJFJF + JJdJFF + JJJJF';
i = 1;
Nall =100;%all possible number of terms
terms = cell(Nall,1);
mults = ones(Nall,1);
Nterms = 0;
buf = '';
flag = 0;
while i <= length(Fstr) + 1
    %parse string
    if i <= length(Fstr)
        if Fstr(i) ~= ' '
            buf = [buf, Fstr(i)];
        else
            flag = 1;
        end
    else
        flag = 1;
    end
    if flag
        %buf is complete, parse it
        buf2 = '';
        multstr = '';
        multdigs = 0;
        for k = 1:length(buf)
            if isstrprop(buf(k),'digit')
                multstr = [multstr,buf(k)];
                multdigs = multdigs + 1;
            end
            if buf(k) == 'F'
                buf2 = [buf(multdigs + 1:k-1),'JF',buf(k+1:end)];
                terms{Nterms + 1,1} = buf2; %add buf to string
                if multdigs > 0
                    mults(Nterms + 1) = str2double(multstr);
                end
                Nterms = Nterms + 1;
                
            end
            if buf(k) == 'J'
                buf2 = [buf(multdigs + 1:k-1),'dJF',buf(k+1:end)];
                terms{Nterms + 1,1} = buf2; %add buf to string
                if multdigs > 0
                    mults(Nterms + 1) = str2double(multstr);
                end
                Nterms = Nterms + 1;
            end
        end
        
        buf = '';
        i = i + 2; %skip "+ "
        flag = 0;
    end
    i = i + 1;
end

skip = zeros(1,Nterms);

dFstr = '';
dFstrfun = '';

for i = 1:Nterms
    if skip(i) == 0
        %find similar
        ctr = mults(i);
        
        for j = i + 1:Nterms
            if strcmp(terms{i,1},terms{j,1})
                ctr = ctr + mults(j);
                skip(j) = 1;
            end
        end
        if ctr > 1
            dFstr = [dFstr, num2str(ctr)];
            dFstrfun = [dFstrfun, num2str(ctr),'*'];
        end
        tmp = terms{i,1};
        dFstr = [dFstr, tmp];
        dFstrfun = [dFstrfun, tmp(1:end-1)];
        if i < Nterms
            dFstr = [dFstr,' + '];
            dFstrfun = [dFstrfun,' + '];
        end
    end
end
dFstr
dFstrfun