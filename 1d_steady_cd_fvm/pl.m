fileID = fopen('result.txt');
step = 5;
n = 102;
y = zeros(n);
x = zeros(n);
for i = 2:(n-1)
	x(i) = 0.005 + 0.01*(i-2);
end
x(1) = 0;
x(102) = 1;
y(1) = 100;
y(102) = 0;
for i = 1:step
	for j =2:101
		y(j) = fscanf(fileID,'%f',1); 
	end
	plot(x,y,"linewidth",2,"r");
	hold on;						
end
