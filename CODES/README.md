# Map-Based-Colocation-Mining
This github includes the code for Abigail Kelly's Qualifying Exam.  The following files are included:
* regional_colocation_compression: this folder contains the code for our map-based regional colocation mining framework.
   * c_functions.cpp: c++ functions used in the regional colocation framework
   * distance_threshold_calc.ipynb: python code that estimates the optimal spatial neighborhood relationship constraint
   * regional_colocation.ipynb: python code that calls the c++ code in c_functions.cpp to perform the colocation mining
   * required_files: this folder holds the intermediate data produced by the framework
   * real_world_data: this folder contains 3 real-world data sets along with their shapefiles.  It is recommended that you use the NorthAmerica data set if you are to run the code due to its shorter run time.
* colocation_compression: this folder contains the code for our map-based colocation code (same code as above minus the regional part- this code is only useful for testing sythetic data sets)
   * c_functions.cpp
   * compression.ipynb
   * synthetic data: this folder contains 5 synthetic data sets varying in clumpiness.  They are titled TestCase1_#.csv where # represents the clumpiness.  For example, TestCase1_1.csv is clumpiness 1.


### How to Configure
1. Open the files in the required_files folder and ensure that all the folders are empty
2. Open distance_threshold.ipynb
3. If running the code using the provided data set, ensure that the 1st line in the 2nd cell is set to: directory = 'data'.  If using a different data set, set the variable named directory equal to the respective path.
4. Open regional_colocation.ipynb
5. The 6th line in the 2nd cell is a user-defined prevalence threshold.  This variable can be changed to include more or less prevalent patterns
6. If running the code using the provided data set, ensure that the 8th line in the 2nd cell is set to: shapefile_path = 'data/shapefile' and that the 9th line in the 2nd cell is set to: directory_path = 'data'.  If using a different data set, set the variable named shapefile_path to the path of the shapefile and the variable named directory_path to the path of the data set

### How to Compile and Run
1. Change your current directory to the directory containing c_functions.cpp
2. Open distance_threshold.ipynb and run all the cells
3. Run the following command in the terminal: **g++ -O3 -shared -o c_functions.so -fPIC c_functions.cpp**
4. Open regional_colocation.ipynb and run all the cells
