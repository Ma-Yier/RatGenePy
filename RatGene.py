import scipy.io as sco
import os
import argparse
import numpy as np


def getBiomassID(model.rxns):
    return None


if __name__ == "__main__":
    # load .mat data
    dataFolder = os.path.join(os.getcwd(), "data")
    modelNames = os.listdir(os.path.join(dataFolder))
    modelChoice = [f.split(".")[0] for f in modelNames]

    # parse params
    parser = argparse.ArgumentParser(description='Hyper Parameters for RatGene')
    parser.add_argument(
        'model', 
        type=str, 
        help='name of the .mat model file downloaded from BiGG, should be inside of ./data/',
        choices=modelChoice,
        )
    parser.add_argument(
        'targeMet', 
        type=str, 
        help='The target metabolite, should be a string in mets not metNames',
        )
    parser.add_argument(
        '--biomass', 
        type=str, 
        help='The biomass reaction in the model, default the growth reaction in the model',
        default=None,
        )
    parser.add_argument(
        '--carbon', 
        type=str, 
        help='The carbon source of the defined problem, default EX_glc__D_e',
        default='EX_glc__D_e',
        )
    parser.add_argument(
        '--oxygen', 
        type=str, 
        help='The input oxygen of the defined problem, default EX_o2_e',
        default='EX_o2_e',
        )
    parser.add_argument(
        '--LBbiomass', 
        type=float, 
        help='The lower thresold of the biomass reaction, default 0.05',
        default=0.05,
        )
    parser.add_argument(
        '--LBcarbon', 
        type=float, 
        help='The lower threshold of the input carbon source exchange reaction, default -15',
        default=-15,
        )
    parser.add_argument(
        '--LBoxygen', 
        type=float, 
        help='The lower threshold of the input oxygen source exchange reaction, default -15',
        default=-15,
        )
    parser.add_argument(
        '--maxLoop', 
        type=int, 
        help='The maximum iterations assigned to RatGene ratio-based procedure, default 1000',
        default=1000,
        )
    parser.add_argument(
        '-t', '--timeLimit', 
        type=int, 
        help='The maximum computation time for the method, default NO TIME LIMIT',
        default=float('inf'),
        )
    parser.add_argument(
        '-p', '--pool', 
        type=int, 
        help='The number of solutions obtained from the IBM CPLEX solution pool, default 10',
        default=10,
        )
    parser.add_argument(
        '--type', 
        type=str, 
        help='The type of modification strategy, gene or reaction, default gene',
        default='gene',
        )
    parser.add_argument(
        '-s', '--size', 
        type=bool, 
        help='Whether reduce the size of the output strategy, True or False, default True',
        default=True,
        )
    args = parser.parse_args()
    
    # form model
    modelName = args.model
    readm = sco.loadmat(''.join(modelName,'.mat'), struct_as_record=False) # read Matlab struct as dict
    model = readm[modelName][0][0]
    
    # get target ID
    target = args.targeMet
    try:
        targetName = np.where(model.mets == target)[0][0]
    except IndexError:
        raise NameError("Check the target name, should be inside of model.mets")
    targetID = np.where(model.rxns == "EX_"+target+"_e")[0][0]
    
    # get biomass ID
    biomassName = args.biomass
    if not biomassName:
        biomassID = getBiomassID(model.rxns)
    else:
        try:
            biomassID = np.where(model.rxns == biomassName)[0][0]
        except IndexError:
            raise NameError("Check the biomass name, should be inside of model.rxns")
        
    # get carbon ID and oxygen ID
    carbon = args.carbon
    try:
        carbonID = np.where(model.rxns == carbon)[0][0]
    except IndexError:
        raise NameError("Check the carbon name, should be an exchange reaction inside of model.rxns")
    oxygen = args.oxygen
    try:
        oxygenID = np.where(model.rxns == oxygen)[0][0]
    except IndexError:
        raise NameError("Check the oxygen name, should be an exchange reaction inside of model.rxns")
    
    # assigen lbs to biomass carbon oxygen
    model.lb[biomassID][0] = args.LBbiomass
    model.lb[carbonID][0] = args.LBcarbon
    model.lb[oxygenID][0] = args.LBoxygen
    
    print("suspend")



