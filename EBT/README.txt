Recommend using VM to rebuild the following system environment:
Win7 (32), in order to run the binary code.

1. To run the tracker "EBT.exe", verify the MATLAB Compiler Runtime (MCR) is installed and ensure you have installed version 8.0 (R2012b). It is available in:
http://www.mathworks.com/products/compiler/mcr/index.html 

1) Add the runtime dlls path into the system "PATH"
variable (otherwise the tracker will report dlls not found):
"MATLAB Compiler Runtime\v80\runtime\win32" ==> "PATH"
2) Make sure "MATLAB Compiler Runtime\v80\bin\win32" is in the "PATH" as well.

2. C++ Redistributable Package 2010 is needed and should be available in:
http://www.microsoft.com/download/en/details.aspx?id=26999

3. Externel libraries: Eigen3, OpenCV 2.31. 
(dlls included)

4. -libEdgeBoxGenerateProposals.dll
   -libEdgeBoxGenerateProposals.h
   -libEdgeBoxGenerateProposals.lib     
They (included) contain some matlab functions obtained using the MATLAB compliler.

5. This implementation is not able to use TRAX and it has to manually disable the TRAX test stage, since the testing just freezes there. 
Put the EXE path into the tracker configuration m file. 
No argument is needed.

6. The tracker sometimes crashes for no reason. It appears pretty random.

7. If there is any problem, please contact: gao.zhu@anu.edu.au.






