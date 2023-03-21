function f = move(a,b,c)
temp = a;
if abs(b) < length(a)
    if strcmp(c,'replace') == 1
        if b>0
            temp(b+1:length(a)) = a(1:length(a)-b);
            temp(1:b) = a(length(a)-b+1:length(a));
        end
        if b<0
            b = abs(b);
            temp(1:length(a)-b) = a(b+1:length(a));
            temp(length(a)-b+1:length(a)) = a(1:b);
        end
    end

    if strcmp(c,'move') == 1
        if b>0
            temp(b+1:length(a)) = a(1:length(a)-b);
            temp(1:b) = 0;
        end
        if b<0
            b = abs(b);
            temp(1:length(a)-b) = a(b+1:length(a));
            temp(length(a)-b+1:length(a)) = 0;
        end
    end
    
    if strcmp(c,'tail') == 1
        if b>0
            temp(b+1:length(a)) = a(1:length(a)-b);
        end
        if b<0
            b = abs(b);
            temp(1:length(a)-b) = a(b+1:length(a));
        end

    end
    
    if strcmp(c,'shuffle') == 1
        temp = shuffle(a);
        if b>0
            temp(b+1:length(a)) = a(1:length(a)-b);
        end
        if b<0
            b = abs(b);
            temp(1:length(a)-b) = a(b+1:length(a));
        end

    end
    
    
end
if abs(b) >= length(a)
    temp(1:length(temp)) = 0;
end


f = temp;



