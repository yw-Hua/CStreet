# CStreet: a computed *C*ell *S*tates *tr*ajectory inf*e*r*e*nce method for *t*ime-series single-cell RNA-seq data 

### Overview
CStreet is a cell states trajectory construction method for time-series single-cell RNA-seq data. It is written in python (python 3.6 or higher) and is available as a commend line tool and a python library to meet the needs of different users.

[Figure1]

CStreet takes advantage of time-series information to construct the connections of *k*-nearest neighbors within and between time points. Then CStreet calculated the connection probabilities of cell states and visualized the trajectory which may include multiple starting points and paths using a force-directed layout method. 

### Installation
CStreet has been packaged and uploaded to PyPI. CStreet and its relevant packages can be installed using one single commands as follows.
   ```shell
   $ pip3 install cstreet 
      # pip3 is the package installer for Python. 
      # If you don't have pip3 on your machine, try install it [https://pip.pypa.io/en/stable/].
   ```
Type the following command to check whether CStreet has been installed successfully.
   ```shell
   $ CStreet -h
   ```

### Quick Start
**Input**: 
   - Expression data: Expression matrix containing the time-series expression level as reads counts or normalized values in tab delimited format, and anndata format are accepted as the input of CStreet. (For example: ExpressionMatrix_t1.txt,ExpressionMatrix_t2.txt,ExpressionMatrix_t3.txt)
   - Cell states info: The cell states information can be inputted by the user or generated using the internal clustering function of CStreet. (For example: CellStates_t1.txt,CellStates_t2.txt,CellStates_t3.txt)

**Command line tool**
   ```shell
   $ CStreet -i ExpressionMatrix_t1.txt,ExpressionMatrix_t2.txt,ExpressionMatrix_t3.txt -s CellStates_t1.txt,CellStates_t2.txt,CellStates_t3.txt -n YourProjectName
   ```
   
**Output**: 
   - An visulization of inferred cell states trajectory. (YourProjectName.CellStatesTrajectory.pdf)
   - The clustered cell states information if not provided by users. (YourProjectName.CellStates.txt)
   - The connection probabilities of cell states. (YourProjectName.ConnectionProbabilities.txt)

[SFig1]

**An example of inferenced cell trajectory**:

![results.png](https://github.com/yw-Hua/MarkdownPicture/blob/master/CStreet/results2.png?raw=true)

