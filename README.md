# SEquence CApture PRocessor (SECAPR)

![SECAPR](https://raw.githubusercontent.com/AntonelliLab/seqcap_processor/master/images/secapr_logo_small.png)

[![downloads](https://anaconda.org/bioconda/secapr/badges/downloads.svg)](http://bioconda.github.io/recipes/secapr/README.html)

**Original Publication: [https://doi.org/10.7717/peerj.5175](https://doi.org/10.7717/peerj.5175)**

<br/>
<br/>

## **Documentation: [Click here](http://antonellilab.github.io/seqcap_processor/)**

<br/>
<br/>
<br/>

***

## Installation & Setup
SECAPR is available as a conda package on the bioconda channel. This makes installation of the SECAPR pipeline with all dependencies very simple. Follow the instructions on this page to get the SECAPR pipeline set up and ready to use.

<br/>

***

### 1. Set up conda

 Conda is a software and environment manager, that makes installation of new software and of required dependencies very simple and straightforward. If you don't already have conda installed on your computer, you can rather quickly and easily install the light-weight version **Miniconda**.

Download the **Python3 version** of [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for your operating system.

**Windows** users can download the `exe` file and install it by double clicking on the file and following the instructions.

If you are a **Mac** user, choose the `pkg` conda download file for installing conda using a graphic user interface. Alternatively you can download the `bash` file (same for **Linux**). In that case you need to open a command line window, navigate to the folder containing the miniconda download, and install it with `sh Miniconda3-*.sh`.

Once installed you can use conda via a bash command line terminal. **Mac** and **Linux** users can simply use the default Terminal software. **Windows** users are recommended to use the **Anaconda Powershell Prompt**, which is automatically installed with Miniconda, instead of the regular Command Prompt.

From your bash-command line terminal add the following channels to your conda installation manager. These channels are the locations where conda finds all required software dependencies.

```bash
conda config --add channels defaults;
conda config --add channels conda-forge;
conda config --add channels bioconda;    
```


***

### 2. Install the SECAPR environment

Besides enabling easy installation, conda also functions as an extremely light-weight **virtual environment** manager. This means you create a virtual environment on your system and install software within that environment without interferring with your main system. For example if you download a conda package that has a specific R version as dependency, you can download it in a virtual environment and it will not affect your current R version that you are using on computer. You can easily connect and disconnect from these virtual environments you create with conda, see below.

Since the SECAPR pipeline has many third-party dependencies, which may effect existing installations of those softwares on your computer, we strongly recommend to install SECAPR in a virtual conda environment. Here we will name this environment `secapr_env`, but you can give it any name you want.

This is super easy with conda. All you have to do to install SECAPR and all its dependencies in a virtual environment is to run the following command:

```bash
conda create -n secapr_env secapr
```

More information about installing SECAPR, including installation of previous versions and using docker containers, can be found by clicking on the badge below:

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/secapr/README.html)

***

### 3. Activate the environment
To activate the newly created `secapr_env` environment, type:

```bash
conda activate secapr_env
```

When the environment is activated, all the necessary software dependencies will be available in the standarad path, e.g. when you type `samtools` the samtools version required for SECAPR will be executed. After you are done using SECAPR, you can deactivate the environment to switch back to your standard environment with this command:

```bash
conda deactivate
```

***

### 4. Check active environment
Check if you are connected to the correct environment (there should eb a star in front of secapr_env in the output of this command):

```bash
conda info --envs
```

In case you are not connected remember to always connect to the `secapr_env` before using SECAPR, by typing:

```bash
conda activate secapr_env
```

To test if the SECAPR installation is working, run the basic help function:

```bash
secapr -h
```

You should see an output with all available SECAPR functions. You can check run the help command for individual functions by typing `secapr name_of_function -h`, e.g.:

```bash
secapr clean_reads -h
```

***

### 5. Install SECAPR development version

The development version of SECAPR is stored on this GitHub page and contains the newest updates, which might not yet be available through the conda version. However you need to install the SECAPR environment with conda first by following the steps above. Once the environment is installed, you can update SECAPR to the development version by following these steps:

1. Connect to your SECAPR environment:

   ```bash
   conda activate secapr_env
   ```

2. [Download](https://github.com/AntonelliLab/seqcap_processor/archive/master.zip) the latest SECAPR version from GitHub  
3. Unzip the downloaded file
4. Move the unzipped directory to a safe location on your computer, i.e. not on your Desktop or Download folder. This will be the path where SECAPR will be executed from in the future.
5. Enter the unzipped SECAPR directory from your command line, replace `/path/to` with the path where you storred the folder:

   ```bash
    cd /path/to/seqcap_processor-master
    ```

6. Remove the installed conda SECAPR version:
   ```bash
   conda remove --force secapr
   ```

7. Install SECAPR from the folder 

    ```bash
    python setup.py install
    ```