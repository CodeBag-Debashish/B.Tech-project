clear all;
clc;
fileID = fopen('temp.out');
n = 100;
res = zeros(100,100);
surf(res);
pause(5)
for step = 1:28
	for col = 1:n
		for row = 1:n
			res(row,col) = fscanf(fileID,'%f',1);
		end	
	end
	surf(res);
	zlim([-0.1 3000])
	xlabel('X');
	xlabel('Y');

	zlabel('Temparature');
	title(['\fontsize{20} Time Step = ' num2str(step)],'Color', 'r','fontweight','bold');
	pause(0.001);
end
