# sonic-sculpting
CS476a Final Project

Project page at https://ccrma.stanford.edu/~mmara/256a/sonic-sculpting/index.md.html

Requires the latest svn version of the G3D innovation engine, available here:

svn://graphics-svn.cs.williams.edu/g3d
user: g3d
pass: g3d

G3D is currently finding a new home after Source Forge, I will update this page when a more permanent home is found.

To compile this project use the provided visual studio 2015 files on windows, or use compile which comes with G3D on OS X (icompile --opt --run)

**User Manual**

- Press the space bar to begin creating a "sculpture" of sound, release the space bar to finish your creation. Your sound will take flight and orbit together with any other sculpture you have created. 
- Press the number keys, '1' meaning the first sculpture you created, and so on, to play back the corresponding sculpture. 
- Press 'Enter' to do something very interesting: play the sounds by sweeping a plane, acting as a playhead, across all the sculptures you see. This reveals the main conceit of the program: the sounds are embedded in 3-dimensional space, and can be resampled using geometric operations.
- Freeze/unfreeze movement with 'f'.
- Ctrl-S to save your soundscape
- Ctrl-L to load it back in
- Press F2 to toggle between free-flight and "on-rails" camera modes.
- , and . make the camera zoon in and out towards the center in "on-rails" camera mode.
