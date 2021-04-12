function [Output] = calculate_5variables(DataInput,option)
% Calculate mean, MSSD1, MSSD2, SD, VSD
    itmp.array1 = DataInput(1:end-1);
    itmp.array2 = DataInput(1+1:end);
    MEAN = mean(DataInput);
    MSSD1 = sqrt((sum((itmp.array2-itmp.array1).^2)/(length(itmp.array1)))).*1000;
    SD = std(DataInput);
    
if option == 1 % In case of raw values
    MSSD2 = sqrt( ((sum((itmp.array2-itmp.array1).^2))/(length(DataInput)-1)) / (sum(DataInput)/length(DataInput)) ).*1000;
    VSD = ( std(abs(itmp.array2-itmp.array1)) / (sum(DataInput)/length(DataInput)) ).*1000;
elseif option == 2 % In case of dynamic connectivity values (including minus values)
    Data2 = DataInput + 1;
    itmp.array1 = Data2(1:end-1);
    itmp.array2 = Data2(1+1:end);    
    MSSD2 = sqrt( ((sum((itmp.array2-itmp.array1).^2))/(length(Data2)-1)) / (sum(Data2)/length(Data2)) ).*1000;
    VSD = ( std(abs(itmp.array2-itmp.array1)) / (sum(Data2)/length(Data2)) ).*1000;
end
    Output = [MEAN, MSSD1, MSSD2, SD, VSD];

end

