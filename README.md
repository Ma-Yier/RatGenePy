# RatGenePy


[RatGene](https://github.com/Ma-Yier/RatGene) Python Version: Growth to production ratio-based design algorithms for constraint-based metabolic networks for growth-coupled production



## Introduction
This is the python version of [RatGene](https://github.com/Ma-Yier/RatGene).

Cite our paper:

+ Ma Y, Tamura T. RatGene: Gene deletion-addition algorithms using growth to production ratio for growth-coupled production in constraint-based metabolic networks[J]. IEEE Transactions on Computational Biology and Bioinformatics, 2025.



## Dependencies

+ Python >= 3.10

+ Numpy

+ Scipy

+ pyscipopt


## RatGenePy

Grammar:

```
python main.py model modelname targetMet targetname [--biomass] [--carbon] [--oxygen] [--LBbiomass] [--LBcarbon] [--LBoxygen] [--maxLoop] [-t --timeLimit] [-p --pool] [--type] [-s --size]
```

Parameters:
+ `model`: The name of the .mat model file downloaded from BiGG, should be inside of ./data/.

+ `targeMet`: The target metabolite, should be a string in mets not metNames.

+ `--biomass`: The biomass reaction in the model, default the growth reaction in the model, default is `None`.

+ `--carbon`: The carbon source of the defined problem, default `EX_glc__D_e`.

+ `--oxygen`: The input oxygen of the defined problem, default `EX_o2_e`.

+ `--LBbiomass`: The lower thresold of the biomass reaction, default `0.05`.

+ `--LBcarbon`: The lower threshold of the input carbon source exchange reaction, default `-15`.

+ `--LBoxygen`: The lower threshold of the input oxygen source exchange reaction, default `-15`.

+ `--maxLoop`: The maximum iterations assigned to RatGene ratio-based procedure, default `1000`.

+ `-t --timeLimit`: The maximum computation time for the method, default `inf`.

+ `-p --pool`: The number of solutions obtained from the IBM CPLEX solution pool, default `10`.

+ `--type`: The type of modification strategy, gene or reaction, default `"gene"`.

+ `-s --size`: Whether reduce the size of the output strategy, True or False, default `True`.


Outputs:

+ `xobj`: The exchange reaction rate for the production of the target metabolite under the condition of applying the output modification strategy to the model..

+ `alpha`: The adopted value of $\alpha$.

+ `knockout`: The name list of knockout strategy indicates which genes or reactions to be knocked out.

Running Example

The following will generate the gene deletion strategy for growth-coupled production of `succ_e` in `e_colic_core` model with maximum uptake rates of glucose and oxygen being 10. 

```
x_obj, alpha, knockout = python ./main.py model "e_coli_core" targetMet "succ_e" --LBcarbon -10 --LBoxygen -10
```

