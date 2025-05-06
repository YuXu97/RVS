hp1 = -30:0.5:30;
hp2 = -1:0.03:1;

[X1,X2] = meshgrid(hp1,hp2);

Y = zeros(length(hp2));
% for i = 1:length(hp1)
     for j = 1:length(hp2)
        Y(j) = nglglklhd([hp(1:4);hp2(j);hp(end)], Psi, y, kernel, M, method);
    end
% end

% surf(X1,X2,Y)
% xlabel('a_1'); ylabel('a_2')
plot(hp2,Y)