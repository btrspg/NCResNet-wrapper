# NCResNet: Noncoding RNA Prediction Based on a Deep Resident Network of RNA Sequences
A deep learning-based method to predict noncoding RNAs of RNA sequences
# Dwonload
git clone https://github.com/abcair/NCResNet.git  
# Installation
Create virtual environments and install dependencies

conda create -n NCResNetEnv python=3.6  
conda activate NCResNetEnv (or source activate NCResNetEnv)  
conda install r  
install.packages("LncFinder",repos="https://cloud.r-project.org/")  
pip install numpy  
pip install pandas  
pip install sklearn  
pip install biopython  
pip install tensorflow  
pip install keras  
pip install rpy2==3.0.1  

# Inputs
Fasta format file

# Outputs
CSV format file

# Command Line usage
python test.py -i ./demo.fasta -m NCResNet.h5 -o result.tsv
