/*
 * Macro to create a hyperstack from 3d raw stacks
 */


setBatchMode(true);
arg_in=split(getArgument," ");

inputFolder = arg_in[0];
outputFolder = arg_in[1];
suffix = arg_in[2];//'.raw';
img = arg_in[3];//"8-bit";
w = arg_in[4];//"100";
h = arg_in[5];//"100";
oS = arg_in[6];//"0"
n = arg_in[7];//"283";
g = arg_in[8];//"0";

processFolder(inputFolder);
run("Concatenate...", "all_open title=[Concatenated Stacks] open");
saveAs("Tiff", outputFolder+File.separator+"Concatenated Stacks.tif");

function processFolder(input) {
	setBatchMode(true); 
	list = getFileList(input);
	
	for (i = 0; i < list.length; i++) {
			
		if(endsWith(list[i], suffix))
			openRawImages(input, list[i]);
	}
}

function openRawImages(input, file) {
	run("Raw...", "open="+input+File.separator+file+" image="+img+" width="+w+" height="+h+" offset="+oS+" number="+n+" gap="+g);
}
