clear all;
clc;
fid = fopen('result.txt');
result = zeros(101,101);
no_of_rows = 101;
no_of_col = 101;
for step = 1:55
    for i=1:101
        for j=1:101
            result(i,j) = fscanf(fid,'%f',1);
        end
    end
	surf(result);
	title(['\fontsize{10} timestep = ' num2str(step)],'Color', 'r','fontweight','bold');
	pause(1)
end

	
	
