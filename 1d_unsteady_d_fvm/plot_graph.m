fileID = fopen('result.txt');
step = 513;
y = zeros(101);
x = [0:0.01:1];
for i =1:139
	for j =1:101
		y(j) = fscanf(fileID,'%f',1); 
	end
	plot(x,y);
	grid on;
	hold on;
end