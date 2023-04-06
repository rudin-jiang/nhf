% matlab code

clear
clc

% see "boysfun.hpp" for more detail
% about BoysFunMaxN and BoysFunSwitch
BoysFunMaxN = 64;
BoysFunSwitch = 200.0;

xSmallInterval = [
      0.0,   1.0;
      1.0,  10.0;
     10.0, 100.0;
    100.0, 190.0;
    190.0, 200.0;
];

xLargeInterval = [
    200.0 + 1.0e-5, 201.0;
     201.0,   210.0;
     210.0,   300.0;
     300.0,  1000.0;
    1000.0, 10000.0;
];

xSmallList1 = linspace(xSmallInterval(1,1), xSmallInterval(1,2), 6000);
xSmallList2 = linspace(xSmallInterval(2,1), xSmallInterval(2,2), 6000);
xSmallList3 = linspace(xSmallInterval(3,1), xSmallInterval(3,2), 6000);
xSmallList4 = linspace(xSmallInterval(4,1), xSmallInterval(4,2), 6000);
xSmallList5 = linspace(xSmallInterval(5,1), xSmallInterval(5,2), 6000);

xLargeList1 = linspace(xLargeInterval(1,1), xLargeInterval(1,2), 6000);
xLargeList2 = linspace(xLargeInterval(2,1), xLargeInterval(2,2), 6000);
xLargeList3 = linspace(xLargeInterval(3,1), xLargeInterval(3,2), 6000);
xLargeList4 = linspace(xLargeInterval(4,1), xLargeInterval(4,2), 6000);
xLargeList5 = linspace(xLargeInterval(5,1), xLargeInterval(5,2), 6000);

smallList = [
    xSmallList1, xSmallList2, xSmallList3, xSmallList4, xSmallList5
];

largeList = [
    xLargeList1, xLargeList2, xLargeList3, xLargeList4, xLargeList5
];

xSmallFilename = "boysfun_test_data_small_x.dat";
xLargeFilename = "boysfun_test_data_large_x.dat";

print_test_data(BoysFunMaxN, smallList, xSmallFilename);
print_test_data(BoysFunMaxN, largeList, xLargeFilename);


function print_test_data(nMax, xList, filename)

    fileID = fopen(filename, "w");
    fprintf(fileID, "%d\n", length(xList));
    fprintf(fileID, "\n");

    for i = 1:length(xList)
        fprintf(fileID, "%.15e  ", xList(i));
        if (mod(i, 10) == 0)
            fprintf(fileID, "\n");
        end
    end

    fprintf(fileID, "\n");

    for i = 1:length(xList)
        x = xList(i);
        for n = 0 : nMax
            fprintf(fileID, "%.15e  ", boysfun(n, x));
        end
        fprintf(fileID, "\n");
    end

    fclose(fileID);
end

