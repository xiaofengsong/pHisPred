# pHisPred
A bioinformatic tool for identifying protein histidine phosphorylation sites

# Installation

The following software should be installed in your cluster or computer before running the pHisPred.py.

*         Python (>= 3.7), https://www.python.org/downloads/.
*         The scikit-learn module, http://scikit-learn.org/stable/install.html.

In most use cases the best way to install Python and scikit-learn package on your system is by using Anaconda(https://www.continuum.io), which is an easy-to-install free Python distirbution and includes more than 400 of the most popular Python packages. Anaconda includes installers(https://www.continuum.io/downloads) for Windows, OS X, and Linux.

# Usage for predicting pHis sites
	Command:
		python pHisPred.py -f input.fa -o output -t euka
	Options:
		-h, --help      show this help message and exit
		-f input files, --file=input files
                	        enter proteins in .fasta (.fa) format.
		-o output files, --out=output files
                	        assign your output file.
		-t class type, --type=class type
                	        eukaryotes: euka (default); prokaryotes: proka.
		-w window size, --window=window size
                	        specific the window size used for extracting the peptides around histidine sites.
		-m trained model, --model=trained model
                	        specific your own trained model.


# Usage for building classification model

	Command:
		python train_model.py -p pHis.fa -n non_pHis.fa -o model.tsv -t proka
	Options:
		-h, --help          show this help message and exit
		-p positive samples, --positive=positive samples
                	            enter peptides around pHis sites in .fasta (.fa) format.
		-n negative samples, --negative=negative samples
                	            enter peptides around non-pHis sites in .fasta (.fa) format.
		-o output files, --out=output files
                	            assign your output file.
		-t class type, --type=class type
                	            eukaryotes: euka (default); prokaryotes: proka.


# Examples
	
	Predicting pHis sites with protein sequences or peptiede sequences centered by His site (> 31 amino acids)
		python pHisPred.py -f examples/euka_proteins.fa -o euka_result.txt -t euka
		python pHisPred.py -f examples/proka_proteins.fa -o proka_result.txt -t proka
		python pHisPred.py -f examples/negatve.fa -o euka_result.txt -t euka
		
	Predicting pHis sites with the model built on your own data 
		python train_model.py -p exampels/positive.fa -n examples/negative.fa -o euka_model.tsv -t euka
		python pHisPred.py -f examples/euka_proteins.fa -o euka_result.txt -m euka_model.tsv


# Author

pHisPred is developed by Jian Zhao (zhao_doctor@hotmail.com). For questions and comments, please contact Jian or submit an issue on github.
