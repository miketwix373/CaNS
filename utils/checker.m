% Field Checker code
path = "/Users/miguel.perez-cuadrado/CaNS/run/data/";
file = "initCond.txt";

dataRaw = table2array(readtable(path+file));
snap(:,:,1) = dataRaw(1:10,1:50);
snap(:,:,2) = dataRaw(11:20,1:50);
snap(:,:,3) = dataRaw(21:30,1:50);

figure()
subplot(1,3,1)
imagesc(snap(:,:,1))
grid minor
colorbar

subplot(1,3,2)
imagesc(snap(:,:,2))
grid minor
colorbar

subplot(1,3,3)
imagesc(snap(:,:,3))
grid minor
colorbar

title('Post - correction')



