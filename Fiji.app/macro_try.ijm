source_dir = getDirectory("Y:\mDrives\storage4\Guillermo\guillermo_functions\trackIngressions\expk0271_wholeEmbryoTest_II\ingressionTracks_i_thr_01\ingressionTest_36586\crops");
stack1=source_dir+"dir1"+File.separator;
        list = getFileList(stack1);
        fileone=stack1 + list[0];
        nummer=list.length;
        run("Image Sequence...", "open=fileone number=nummer starting=1 increment=1 scale=50 file=[] or=[] convert sort"); 