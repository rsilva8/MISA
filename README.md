# MISA: Multidataset Independent Subspace Analysis
This repository contains a MATLAB implementation for Multidataset Independent Subspace Analysis (MISA), using the Kotz distribution.

![MDM](./MDMmodel_new.png#raw=true "General MDM problem")

# Usage

### Step 1: Clone to MISA folder and start Matlab
```
cd ~
git clone git@github.com:rsilva8/MISA.git MISA
cd MISA
matlab
```

Optional: get example datasets (about 0.5GB)
```
cd ~
git clone git@github.com:rsilva8/MISA-data.git MISA-data
cd MISA-data
unzip MISA-data.zip
cd ../MISA
```

### Step 2: Load data for analysis
```
% Add folders to path
addpath('./scripts')
addpath('./scripts/toy_example')

%% Load the true sources (Y) and true mixing matrices (A) to generate the mixtures

% Load sources
load(fullfile('.','scripts','toy_example','jointsourcesMISA_case6.mat'), 'Sgt')
Y = Sgt;

% Load Tetris-String features
load(fullfile('.','scripts','toy_example','A_TetrisString_orthofeat_r001.mat'),'A');

% Generate mixtures
X{1} = A{1}*Y{1};
X{2} = A{2}*Y{2};
```

### Step 3: Define subspace structure
```
% Define the number of datasets (here, the number of modalities)
M = 1:length(X);

S = cell(1,2);           % Cell array: each cell contains a matrix K x C(m).
                         % Each k-th row has 0's and 1's to indicate what
                         % source go within the k-th subspace in dataset m

% Modality 1 = dataset 1
S{1} = [1 0 0 0;... % source 1 into subspace 1
        0 1 0 0;... % source 2 into subspace 2
        0 0 1 1;... % sources 3 and 4 into subspace 3
        0 0 0 0];   % no sources from modality 1 into subspace 4

% Modality 2 = dataset 2
S{2} = [1 0 0 0 0 0;... % source 1 into subspace 1
        0 1 1 0 0 0;... % sources 2 and 3 into subspace 2
        0 0 0 1 0 0;... % source 4 into subspace 3
        0 0 0 0 1 1];   % sources 5 and 6 into into subspace 4 
```

### Step 4: Create MISA object
```
% Define MISA parameters

get_MISA_parameters

% Initialize MISA object

data1 = MISAKRE(w0, M, S, X, ...
                beta, eta, lambda, ...
                gradtype, sc, preX, ...
                REtype, REapproach, RElambda, ...
                REref, REreftype, REreflambda, rC);
```

### Step 5: Run MISA
```
% Run MISA: PRE + LBFGS-B + Nonlinear Constraint + Combinatorial Optimization

execute_full_optimization

% NOTE: toggle lines 41-43 in @utils/getop.m to see detailed optimization iterations (considerably slower)
```

### Step 6: Check results and visualize output
```
% Check results

fprintf("\nFinal MISI: %.4f\n\n", data1.MISI(A))
% typically, a number < 0.1 indicates successful recovery of the sources
```

> Final MISI: 0.0303
```
%% Visualize recovered (mixing) patterns
view_results
```

- Final MISA estimates

![tetris_est](./scripts/toy_example/tetris_estimate.png#raw=true "MISA estimate - Tetris")
![string_est](./scripts/toy_example/string_estimate.png#raw=true "MISA estimate - String")

This demonstrates that *separation occurs up to subspace identification*.
Features retain their dependence *within subspaces*.

# Citation 
- Link to arxiv

# Dependences
- Some examples assume [SPM] toolbox is readily available.
- Some examples assume [GIFT] toolbox is also available.

[SPM]: https://www.fil.ion.ucl.ac.uk/spm/
[GIFT]: http://trendscenter.org/software/gift/

# Folders

## `@MISAK`
- This is the base MISA class.
Use it to instantiate a MISA object which will contain the input data and the methods necessary for source estimation, including combinatorial optimization.

## `@MISAKRE`
- This class inherits from `@MISAK`.
Use it to gain access to Reconstruction Error (RE) functionality.

## `@utils`
- This class implements various useful utility functions.

## `@gsd`
- A class for simulated data generation by simulating sources (Y) and mixing matrices (A).
Also supports user-supplied mixing matrices (A).

## `@gsm`
- A class for simulated mixture generation.
Supports user-supplied sources (Y) and mixing matrices (A).

## `scripts`
- Scripts for running fully simulated and hybrid examples.

## `other_methods`
- Alternate methods evaluated for comparison.

## `results`
- Some results for the examples include in the `scripts` folder.

# Files
- `optprob4.mat`: contains an optimization problem object for use with MATLAB's optimizaion toolbox.
- `validateFirstDerivatives_.m`: slight modification of MATLAB's function to check the accuracy of user-defined derivatives wrt finite differencing.

# Written For
- MATLAB 2017a

# License