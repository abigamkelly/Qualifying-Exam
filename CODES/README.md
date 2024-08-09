# Map-Based-Colocation-Mining
This github includes the code for Abigail Kelly's Qualifying Exam.  The following files are included:
* regional_colocation_compression: this folder contains the code for our map-based regional colocation mining framework.
   * c_functions.cpp: c++ functions used in the regional colocation framework
   * distance_threshold_calcalculation.ipynb: python code that estimates the optimal spatial neighborhood relationship constraint
   * regional_colocation.ipynb: python code that calls the c++ code in c_functions.cpp to perform the colocation mining
   * required_files: this folder holds the intermediate data produced by the framework
   * real_world_data: this folder contains 3 real-world data sets along with their shapefiles.  It is recommended that you use the NorthAmerica data set if you are to run the code due to its shorter run time.
* colocation_compression: this folder contains the code for our map-based colocation code (same code as above minus the regional part- this code is only useful for testing synthetic data sets)
   * c_functions.cpp: c++ functions called in compression.ipynb
   * compression.ipynb: python code that calls the c++ code in c_functions.cpp to perform the map-based colocation mining
   * required_files: intermediate data produced by the code
   * synthetic data: this folder contains 5 synthetic data sets varying in clumpiness.  They are titled TestCase1_#.csv where # represents the clumpiness.  For example, TestCase1_1.csv is clumpiness 1.


### How to Configure
#### regional_colocation_compression
1. Open the files in the required_files folder and ensure that all the folders are empty
2. Open distance_threshold_calculation.ipynb
3. If running the code using the NorthAmerica data set, ensure that the 1st line in the 2nd cell is set to: directory = 'real_world_data/NorthAmerica'.  If using a different data set, change the data set name.  For example, if using the MiddleEast data set, directory = 'real_world_data/MiddleEast'.
4. Open regional_colocation.ipynb
5. The 6th line in the 2nd cell is a user-defined prevalence threshold (prevalence_threshold).  This variable can be changed to include more or less prevalent patterns.
6. If running the code using the NorthAmerica data set, ensure that the 8th line in the 2nd cell is set to: shapefile_path = 'real_world_data/NorthAmerica/shapefile' and that the 10th line in the 2nd cell is set to: directory_path = 'real_world_data/NorthAmerica'.

#### colocation_compression
1. Open the files in the required_files folder and ensure that all the folders are empty
2. Open compression.ipynb
3. The 2nd cell holds the user-defined variables (distance_threshold and prevalence_threshold).  These can be changed accordingly, but if running the code using any of the synthetic data sets, use distance_threshold = 10 and prevalence_threshold = 0.5.
4. If running the code using the TestCase1_1.csv synthetic data set, ensure that the file path in 2nd line in the 3rd cell is set to: "synthetic_data/TestCase1_1.csv".  If using a different synthetic data set, change the data set name.  For example, if using TestCase1_5.csv, the file path would be "synthetic_data/TestCase1_5.csv".

### How to Compile and Run
#### regional_colocation_compression
1. Change your current directory to the directory containing c_functions.cpp
2. Open distance_threshold_calculation.ipynb and run all the cells
3. Run the following command in the terminal: **g++ -O3 -shared -o c_functions.so -fPIC c_functions.cpp**
4. Open regional_colocation.ipynb and run all the cells

#### colocation_compression
1. Change your current directory to the directory containing c_functions.cpp
2. Run the following command in the terminal: **g++ -O3 -shared -o c_functions.so -fPIC c_functions.cpp**
3. Open compression.ipynb and run all the cells
