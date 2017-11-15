# A model for the avascular phase of tumoral growth

A work from 06/2016 by Charles Gaydon and Hélio Wang.

**TL;TR**: during my third year at the INSA Lyon Biosciences Department, we had the freedom to develop a diffusion-reaction based model in order to address a biological question of our choice. 
Based on [1] and [2] (see *Docs* folder), we implemented a **multiscale model** that aimed to replicate the first phase of the tumoral growth. It is two dimensionna and is based on reaction-diffusion equations and on a protein regulatory network.

**Understanding this step is essential to understand the initial aggressivenes of a tumor!**

### 1. The Avascular Tumor Growth problem


At this point, there is no blood vessels (yet !) to irrigate the fast growing tumorous cells with what they need (oxygene, glucose, growth factors, etc). Hence, only diffusion assures the ability of a cell to survive.  We distinguish three types of cells, namely : Quiescent, Proliferative (i.e. tumoral), and necrotic. Indeed, a tumor comes from a deregulated Quiescent cell that became Proliferative. But as it divides to forms a spheroid of proliferative cells with high metabolism (in particular the glycolytic metabolism), the center of this spheroid becomes less accessible to the environnment nutriments and essential growth factors: diffusion is less effective, and wastes starts to accumulate. It is also said that the acidity of the environment rises, which in turn makes it a deadly place for normal (quiescent) cells to thrive.

**In short**: the cells in the middle of the spheroid necrose, which may be essential to maintain a deadly environnment for normal cells, and accelerating the growth of the tumour.

### 2. How to run the code and produce pretty images
Use :

	git clone https://github.com/CharlesGaydon/Avascular_Tumoral_Growth_Model
	cd Avascular_Tumoral_Growth_Model

You can modify model parameters inside the functions, or run them directly (using python 2.7 and no fancy packages) via :	
	
	python Avascular_Growth.py

or
	
	python Avascular_Growth_with_network.py

The propagation of cells and the concentration of cell types are saved in the *History* folder. 
To generate the images of the simulation in the folder *Figures*, simply run :

	python Generate_figures.py prefix_of_results

where *prefix_of_results* is of form "p0.9size15co1vit0.05Tc3" or "NET_p0.9size15co1vit0.05Tc3".

A subset of what can be produced is already in the folder *Figures*. See in paticular the "p0.1 & gene regulation network.gif" file! Our main results and the further mathematical work made by Helio Wang are not part of this directory.

### 3. Acknowledgment

Our model is essentially based on [1] : 

> "At the cellular scale, our model considers cell growth and proliferation, intercellular adhesion, and necrotic cell death. At the subcellular scale, we include a protein expression regulatory network for the control of cell-cycle arrest." 

For this last part, I implemented a boolean regulation network whose definition can be found in the *Docs* folder. 

> "At the extracellular scale, the model considers diffusion, consumption, and production of nutrients, metabolites, growth promoters, and inhibitors.""

[1] : Jiang et al (2005) - A Multiscale Model for Avascular Tumor Growth
[2] : Gatenbi & Gawlinski (1996) - A Reaction-DiffusionModel of Cancer Invasion
