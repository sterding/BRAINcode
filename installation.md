# Instruction for installing prerequisitions of BRAINcode pipeline

## fastq-mcf

```bash
git clone https://github.com/ExpressionAnalysis/ea-utils.git
cd ea-utils/clipper
PREFIX=$HOME make install
fastq-mcf -h  # test if installation is successful
```
If successful, this should install fastq-mcf in `$HOME/bin` folder. If not in system path, try `export PATH=$HOME/bin:$PATH`

## kpal

```bash
git clone https://github.com/LUMC/kPAL.git
cd kPAL
pip install --user -e .
kpal -h # to test if installation is successful
```
    
## UCSC Kent utility

If not installed in your system already, see details here:
http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=blob;f=src/userApps/README