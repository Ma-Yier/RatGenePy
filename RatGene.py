import scipy.io as sco
import os
import argparse
import numpy as np
import time


def solver(objfunc, **kwargs):
    return None
# return 0 if no result

def introExchange(model,biomassID,input_list:list,targetName):
    if len(input_list) == 1:
        carbonID = input_list[0]
        oxygenID = 0
    elif len(input_list) == 2:
        carbonID = input_list[0]
        oxygenID = input_list[1]
    else:
         raise IndexError("Number of input exchange reactions should be 1 or 2.")
    
    ex_target = "EX_" + targetName + "_e"
    dx_target = "DX_" + targetName + "_e"

    # find targetID
    try:
        targetID = np.where(model.rxns == ex_target)[0][0]
    except:
        targetID = None
    try:
        targetID = np.where(model.rxns == dx_target)[0][0]
    except:
        targetID = None
    # check legality of target met
    if targetID == carbonID or targetID == oxygenID:
        raise NameError("target mets cannot be oxygen or carbon source.")
    
    if not targetID:
        targetID = len(model.rxns)
        metID = np.where(model.mets == targetName)[0][0]
        model.rxns = np.append(model.rxns, [ex_target], axis=0) 
        model.S = np.append(model.S, np.zeros(len(model.mets),1), axis=1)
        model.S[metID][targetID] = -1
        model.grRules = np.append(model.grRules, [''], axis=0)
        model.lb = np.append(model.lb, [-1000], axis=0)
        model.ub = np.append(model.ub, [1000], axis=0)
        model.c = np.append(model.c, [0], axis=0)
    
    # calculate TMPR
    model.c[biomassID][0] = 0
    model.c[targetID][0] = 1
    x = solver(-model.c, Aeq=model.S, beq=model.b, lb=model.lb, ub=model.ub)
    TMPR = x[targetID]
    model.c[biomassID][0] = 1
    model.c[targetID][0] = 0

    return model, targetID, TMPR

def RatMethodRxn(model,biomassID,targetID,maxLoop,gap,timeLimit):
    alpha = 0
    knockout = []
    x_target = 0
    m = model.S.shape[0]
    n = model.S.shape[1]
    ratio = [0]*n
    ratio[biomassID] = -alpha
    ratio[targetID] = 1
    sAeqm = np.append(model.S, ratio, axis=0)
    Aeqm = np.concatenate((sAeqm, np.zeros(m+1, n)), axis=1)
    beqm = np.append(model.b, [0], axis=0)
    Am1 = np.concatenate((np.eye(n), -np.eye(n)), axis=1)
    Am2 = np.concatenate((-np.eye(n), -np.eye(n)), axis=1)
    Am = np.concatenate((Am1, Am2), axis=0)
    bm = np.zeros(2*n, 1)
    c = np.concatenate((np.zeros(n, 1), np.ones(n, 1)), axis=0)
    lbm = np.concatenate((model.lb, np.zeros(n, 1)), axis=0)
    ubm = np.concatenate((model.ub, np.ones(n, 1)*999999), axis=0)

    s = time.time()
    for i in range(maxLoop):
        alpha += gap
        Aeqm[m][biomassID] = -alpha
        x = solver(c, A=Am, b=bm, Aeq=Aeqm, beq=beqm, lb=lbm, ub=ubm)
        x_target, knockout = verifyRxn(model, x[:n][:], biomassID, targetID, model.lb[biomassID][0])
        if x_target > 0:
            break
        # check time limit
        e = time.time()
        if e - s > timeLimit:
            break 
    return x_target, knockout


def verifyRxn(model, knot, biomassID, targeID, LBbiomass):
    x_target = 0
    knot = np.abs(knot)
    delist = np.where(knot < 0.0000001)
    model.lb[delist[0]][delist[1]] = 0
    model.ub[delist[0]][delist[1]] = 0
    model.lb[biomassID][0] = 0

    x = solver(-model.c, Aeq=model.S, beq=model.b, lb=model.lb, ub=model.ub)
    if x[biomassID] > LBbiomass:
        x_target = x[targeID]

    return x_target, zip(delist[0], delist[1])


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
    
    # get biomass ID
    biomassName = args.biomass
    if not biomassName:
        biomassID = np.nonzero(model.c)[0][0]
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

    # check target exchange rxns exist or not
    model, targetID, TMPR = introExchange(model,biomassID,[carbonID, oxygenID],targetName)
    if TMPR <= 0:
        raise ValueError("Theoritcal Maximum Production Rate is not positive.")
    # RatMethod for genes or rxns
    gap = TMPR/(args.maxLoop*args.LBbiomass)
    if args.type == "gene":
        RatMethodGene()
    elif args.type == "reaction":
        RatMethodRxn(model,biomassID,targetID,args.maxLoop,gap,args.timeLimit)
    else:
        raise NameError("Deletion type should be either \'gene\' or \'reaction\'")


    print("suspend")



