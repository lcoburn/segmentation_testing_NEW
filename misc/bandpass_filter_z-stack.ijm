
setBatchMode(true);

//o_img=File.openDialog("Select a stack");
arg_in=split(getArgument," ");
o_img=arg_in[0];
f_img=arg_in[1];
open(o_img);

name=getTitle();

rename("file_in")

selectWindow("file_in");

run("Bandpass Filter...", "filter_large=20 filter_small=1 suppress=None tolerance=5 autoscale saturate process");
run("Gaussian Blur...", "sigma=1 stack");


//f_path=File.getParent(f_img);
f_name=File.getName(f_img);

rename(f_name);
//		rename("filtered_"+name);

setBatchMode(false);

//saveAs("Tiff", "/mnt/cjwsmb/Antti/time_sequences_2d_tracking/No122_241111/segmentation_revisited/filtered_MA1_"+name);
//saveAs("Tiff", File.getParent(o_img)+"/filtered_"+name);
saveAs("Tiff", f_img);
