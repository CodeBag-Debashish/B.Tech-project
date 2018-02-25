clear all;
clc;
fileID = fopen('temp.out');
n = 100;
res = zeros(100,100);

for step = 1:1
	for col = 1:n
		for row = 1:n
			res(row,col) = fscanf(fileID,'%f',1);
		end	
	end
	surf(res);
	xlabel('X');
	xlabel('Y');
	zlabel('Temparature');
	

	%zlim([0 500]);
	%hold on;
	%title(['\fontsize{16} generation = '],'Color', 'b','fontweight','bold');
	%if step == 1
	%	pause(5);
	%end
	%pause(2);
end
