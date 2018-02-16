The timbre space example requires GEM for plotting and Pd version 0.47-1 or greater.

A few things have changed since the previous timbre space example patch. First, this version now runs audio as a subprocess using [pd~]. This puts the audio and video components on separate threads to avoid audio dropouts when graphics processing gets intense. You will not see the audio process patch unless you choose to load it without the -nogui option, or open it directly from the 09-timbre-space/lib/ directory.

Second, in order to avoid dependency on [split_path] and [folder_list], this version requires a list of file names specified in a text file. Some example file lists are provided in the filename-lists directory. They list audio file names relative to the 09-timbre-space/lib/ directory, where the audio patch resides. Another option is to simply provide the full path to your file names so that you don't need to copy or move audio files to the 09-timbre-space/audio/ directory. To automatically write full file paths to a text file, cd to a directory containing audio files, and run:

ls -d -1 $PWD/*.wav > filenames.txt

This will produce filenames.txt, where each line shows the full path of each .wav file in the directory. Remember that any such text files containing audio file names must be in the 09-timbre-space/filename-lists/ directory for the patch to find them.

Finally, this version offers a mode where mousing over points can play back entire sound files rather than grains of a fixed window size. Activate the "full-sample-mode" toggle for that, and choose a fade in/out time in milliseconds (20ms by default). The window size and dB threshold settings do not take effect when using full sample playback mode.

After you have loaded and analyzed your sounds, active GEM rendering and the "Mouseover" toggle to enable sound playback when you mouse over any point. Adjust the mouse radius slider to change the mouse's range of influence. Use the zoom and offset controls to navigate to specific parts of the plot and get a more detailed look.

To assign different features to the X Y and Z axes, just click the radio buttons at the top of the patch.  Note that most of the new audio features available in timbreID 0.7 have been added to this patch. Mouse-over audio browsing is only possible in 2 dimensions, but you can also try plotting in 3 dimensions and rotate the space to get a better look at things. If you've got further questions or suggestions, send me an email: w@williambrent.com