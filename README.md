# plumeria_wd
One-dimensional model for wind-driven volcanic plumes, This is the code used for analysis in Mastin (2014).

Reference:
Mastin, L.G., 2014. Testing the accuracy of a 1-D volcanic plume model in estimating mass eruption rate. Journal of Geophysical Research: Atmospheres, 119(5).

This program was written using Linux CentOS 7, but should work on all linux or linux-like operating  systems (e.g. ios, Unix).  

To run the model, you need a fortran compiler.  The free gnu fortran compiler runs on all linux machines and can be downloaded in CentOS using the command:

   sudo yum install gfortran
	
In the makefile, you may need to modify lines 13-15, and 23-36 in the makefile depending on the compiler you are using.
   
INSTALLATION INSTRUCTIONS:
<p>After cloning the repository to your local computer, unzip the repository, compile the program, and try running it using the following  commands (or equivalents for your OS).  Comments explaning these command are listed to the right, and start with a hashtag (#)<br>
    unzip plumeria_wd                       #unzips the zipped file<br>
	 cd plumeria_wd                          #changes directories to the new directory of unzipped files<br>
	 [edit plumeria_wd/makefile]             #Edit makefile to ensure that it is using the appropriate compiler<br>
	 make                                    #compiles the program<br>
	 ./plumeria_wd input/default.inp         #runs the program using the default input values stored in the file input/default.inp.<br>
</p?	 
    
After the test run, the output file "output/default_params_out.txt" should be identical to the file of the same name in the example_output directory.    
