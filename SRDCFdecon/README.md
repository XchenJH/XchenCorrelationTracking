This MATLAB code implements the SRDCFdecon tracker [1], which employs the unified discriminative tracking formulation [1] for decontaminating the training set in the SRDCF [2].

Installation:
If the pre-compiled mexfiles do not work, the provided compilemex_[linux|windows] should compile them from the provided source code. Alternatively, you can try to modify them for your system.

Instructions;
* The "demo.m" script runs the tracker on the provided "Couple" sequence.
* The runfiles/ folder contain default parameter settings for the OTB and TempleColor datasets as used in the paper.
* The the runfiles supports the OTB dataset interface.
* The "test.m" runfile contain the OTB_settings by default.
* Debug visualization can be turned on or off in the runfiles (e.g. in test.m).

Project webpage:
http://www.cvl.isy.liu.se/research/objrec/visualtracking/decontrack/index.html

Contact:
Martin Danelljan
martin.danelljan@liu.se

Third party code used in the implementation of this tracker is:
* Piotrs image processing toolbox [3]
* mtimesx [4]
* opencv [5]
* lightspeed toolbox [6]

[1] Martin Danelljan, Gustav Häger, Fahad Khan, Michael Felsberg. 
	Adaptive Decontamination of the Training Set: A Unified Formulation for Discriminative Visual Tracking.
	In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2016.

[2] Martin Danelljan, Gustav Häger, Fahad Shahbaz Khan and Michael Felsberg.
	Learning Spatially Regularized Correlation Filters for Visual Tracking.
	In Proceedings of the International Conference in Computer Vision (ICCV), 2015. 

[3] Piotr Dollár.
    "Piotr’s Image and Video Matlab Toolbox (PMT)".
    http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html

[4] http://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support

[5] http://opencv.org/

[6] http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
