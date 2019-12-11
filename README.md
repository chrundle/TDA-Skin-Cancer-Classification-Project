# MTG 7396 - Topological Data Analysis Project
### Classification of Skin Cancer Images using Topological Data Analysis

# Github Contents

Here we describe the directories included in the github repository for this project:

* **Code**: Contains files for preprocessing data and training classification models.
* **Documents**: Contains .pdf copy of the project presentation and report and subdirectories with .tex files for compiling the presentation and report.

### Code

This directory contains all of the code that was used to preprocess the original [Skin Cancer MNIST: HAM 10000](https://www.kaggle.com/kmader/skin-cancer-mnist-ham10000) data set. There are three main files that were used to preprocess the data files downloaded from Kaggle. Details of these files are provided below.

1. `Preprocess_Images_Sampling_Routine.ipynb`: Includes routines for sampling points from skin cancer images and saves sampled point data to .csv file so that it can be loaded in **R** for topological data analysis.

2. `DeathVector_PersLand_SimpleSampled.R` and `DeathVector_PersLand_ClusterSampled.R`: Either of these files should be used after the image preprocessing and sampling routine file. These files load the saved .csv data from the image sampling routine and compute the death vectors and persistence landscape vectors for the sampled point clouds. The first file creates the death vectors and persistence landscape vectors for a point cloud generated using the simple sampling method while the second file creates the death vectors and persistence landscape vectors for a point cloud generated using the cluster sampling method. The main difference is the line where the .csv files are loaded. The majority of these files are the same but were kept separate to avoid any confusion or unintentional mixing of the data. These files also have routines implemented for computing a principal component analysis of the resulting death vectors and persistence landscapes and use the ksvm pacakge in **R** to train a SVM classifier for the death vectors/persistence landscape vectors.

The remaining file in the **Code** directory contains the deep neural network classifier and is described below.

3. `NeuralNetworkClassifiers.ipynb`: Includes routine for loading saved death vectors and persistence landscapes from the `*.R` files and converts them to a format for use with TensorFlow. There are various neural network topologies in place for classifiers on the death vectors or on the persistence landscapes.

### Documents

This directory contains two files that summarize the key details from the project:

* `Report.pdf`: A report outlining more details on the project, key takeaways from the project, and future work that could be done in this area.

* `Presentation.pdf`: An in-class presentation that was given for the project.
