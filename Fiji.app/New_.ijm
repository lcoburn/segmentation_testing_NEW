/*
 * Macro to create a hyperstack from 3d raw stacks
 */


inputFolder = "Y:\\mDrives\\storage4\\Guillermo\\guillermo_functions\\trackIngressions\\expk0271_wholeEmbryoTest_II\\ingressionTracks_i_thr_01\\ingressionTest_36586\\crops";
outputFolder =  "Y:\\mDrives\\storage4\\Guillermo\\guillermo_functions\\trackIngressions\\expk0271_wholeEmbryoTest_II\\ingressionTracks_i_thr_01\\ingressionTest_36586\\crops";

suffix = '.raw';
img = "8-bit";
w="100";
h="100";
oS="0"
n="283";
g="0";

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
	// do the processing here by replacing
	// the following two lines by your own code
	print(File.separator);
	run("Raw...", "open="+input+File.separator+file+" image="+img+" width="+w+" height="+h+" offset="+oS+" number="+n+" gap="+g);

}
