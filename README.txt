> In this folder all scripts uses relative paths to call the functions.

> Run "runTracking.m" using F5 to perform segmentation and tracking,
  results in output directory.

> IMPORTANT: for segmentation a bandpass filter is applyied to the images
  using imageJ/Fiji. There is a macro for that in the "misc" directory.
  I have included Fiji binaries in the folder, but maybe it is necessary
  add them to the path.

> EVEN MORE IMPORTANT: To call Fiji from matlab the script "evaluate_imagej_script.m"
  is used. It contains the commands to call Fiji from the command line in Windows.
  For MAC OS or LINUX it MUST BE MODIFIED!!. The function to call Fiji that
  I normally use ("run_imagej_script_windows") is included for inspiration.