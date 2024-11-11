clc; clear;
root = pwd;
path = root + "/run/data/";

filenames = ["bl_init_u.out" "bl_init_v.out" "bl_init_w.out"];

for i=1:3
    data = readtable(path+filenames(i),"FileType","text");
    temp = table2array(data);
    profile(:,i)= temp(:,2);
    z = temp(:,1);
end

eta = z/0.02
figure()
plot(profile(:,1),z)