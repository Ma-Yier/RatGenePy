import scipy.io as sco
import os
import argparse
import numpy as np
import time
from pyscipopt import Model, quicksum
import copy

def solver(objfunc, **kwargs):
    if "opt" in kwargs:
        pass
    else:
        raise KeyError("should assign optimization type LP or MILP")
    
    
    model = Model()
    # var nums and add vartypes
    n = kwargs["Aeq"].shape[1]
    if kwargs["opt"] == "lp":
        vartype = ["C"] * n
    elif kwargs["opt"] == "milp":
        vartype = kwargs["vtype"]
    else:
        raise ValueError("opt should be either lp or milp")

    # add vars
    lb = kwargs["lb"]
    ub = kwargs["ub"]
    x = []
    for i in range(n):
        x.append(model.addVar(name=f"x_{i}", vtype=vartype[i], lb=lb[i], ub=ub[i]))

    # add cons
    y = []
    Aeq = kwargs["Aeq"]
    beq = kwargs["beq"]
    for i in range(Aeq.shape[0]):
        expr = sum(Aeq[i][j] * x[j] for j in range(n))
        y.append(model.addCons(expr == beq[i]))
    
    if "A" in kwargs and "b" in kwargs:
        A = kwargs["A"]
        b = kwargs["b"]
        for i in range(A.shape[0]):
            expr = sum(A[i][j] * x[j] for j in range(n))
            y.append(model.addCons(expr <= b[i]))

    if "pool" in kwargs:
        model.setParam("limits/soluions", kwargs["pool"])

    # obj func
    model.setObjective(quicksum(objfunc[i][0] * x[i] for i in range(n)), sense="minimize")
    model.optimize()

    return model, x
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
        model.rxns = np.append(model.rxns, np.array([[ex_target]]), axis=0) 
        model.S = np.append(model.S, np.zeros((len(model.mets),1)), axis=1)
        model.S[metID][targetID] = -1
        model.grRules = np.append(model.grRules, np.array([['']]), axis=0)
        model.lb = np.append(model.lb, np.array([[-1000]]), axis=0)
        model.ub = np.append(model.ub, np.array([[1000]]), axis=0)
        model.c = np.append(model.c, np.array([[0]]), axis=0)
    
    # calculate TMPR
    model.c[biomassID][0] = 0
    model.c[targetID][0] = 1
    optres, index = solver(-model.c, Aeq=model.S, beq=model.b, lb=model.lb, ub=model.ub, opt="lp")
    x = optres.getSols()[0]
    TMPR = x[index[targetID]]
    model.c[biomassID][0] = 1
    model.c[targetID][0] = 0

    return model, targetID, TMPR

def RatMethodRxn(model,biomassID,targetID,maxLoop,gap,timeLimit):
    alpha = 0
    knockout = []
    x_target = 0
    (m, n) = model.S.shape
    ratio = np.array([[0]*n])
    ratio[0][biomassID] = -alpha
    ratio[0][targetID] = 1
    sAeqm = np.append(model.S, ratio, axis=0)
    Aeqm = np.concatenate((sAeqm, np.zeros((m+1, n))), axis=1)
    beqm = np.append(model.b, np.array([[0]]), axis=0)
    Am1 = np.concatenate((np.eye(n), -np.eye(n)), axis=1)
    Am2 = np.concatenate((-np.eye(n), -np.eye(n)), axis=1)
    Am = np.concatenate((Am1, Am2), axis=0)
    bm = np.zeros((2*n, 1))
    c = np.concatenate((np.zeros((n, 1)), np.ones((n, 1))), axis=0)
    lbm = np.concatenate((model.lb, np.zeros((n, 1))), axis=0)
    ubm = np.concatenate((model.ub, np.ones((n, 1))*999999), axis=0)

    s = time.time()
    for i in range(maxLoop):
        alpha += gap
        Aeqm[m][biomassID] = -alpha
        opt, index = solver(c, A=Am, b=bm, Aeq=Aeqm, beq=beqm, lb=lbm, ub=ubm, opt="lp")
        x = opt.getSols()[0]
        x_target, knockout = verifyRxn(model, [opt.getSolVal(x, value) for value in index[:n]], biomassID, targetID, model.lb[biomassID][0])
        if x_target > 0:
            break
        # check time limit
        e = time.time()
        if e - s > timeLimit:
            break 
    return x_target, alpha, knockout


def cmpMatDiff(filename, matrix):
    mat_data = sco.loadmat("cache/"+filename+".mat")
    mat_matrix = mat_data[filename]

    same_shape = matrix.shape == mat_matrix.shape
    same_values = np.allclose(matrix, mat_matrix, atol=1e-8)

    print("same shape:", same_shape)
    print("same value:", same_values)

    if same_shape and same_values:
        print("two matrix equal")
    else:
        print("different")

    diff_mask = np.abs(matrix - mat_matrix) > 1e-8
    diff_indices = np.argwhere(diff_mask)
    print("following entry different:")
    for idx in diff_indices:
        idx_tuple = tuple(idx)
        print(f"position {idx_tuple}: a={matrix[idx_tuple]}, b={mat_matrix[idx_tuple]}")
    return diff_indices
    

def RatMethodGene(model,biomassID,targetID,TMGR,maxLoop,gap,numMultiSols,timeLimit):
    lessMatrixLeft, lessMatrixRight, equalMatrixLeft, equalMatrixRight, nGpr, nAux, nGen, _ = constructMatrix(model)
    alpha = 0
    n = model.S.shape[1]
    row_ratio = equalMatrixLeft.shape[0]

    # construct A and b
    Am = lessMatrixLeft
    bm = lessMatrixRight

    # construct Aeq and beq
    ratioCons = np.zeros((1,n+nGpr+nAux+nGen))
    ratioCons[0, targetID] = 1
    ratioCons[0, biomassID] = 0
    Aeqm = np.concatenate((equalMatrixLeft, ratioCons), axis=0)
    beqm = np.concatenate((equalMatrixRight, np.array([[0]])), axis=0)

    # set lb and ub
    lb = np.concatenate((model.lb, np.zeros((nGpr+nAux+nGen,1))), axis=0)
    ub = np.concatenate((model.ub, np.ones((nGpr+nAux+nGen,1))), axis=0)

    # set object function 
    fobj = np.concatenate((np.zeros((n,1)), TMGR*np.ones((nGpr,1)), np.zeros((nAux+nGen,1))), axis=0)
    fobj[biomassID, 0] = -1

    # set milp labels
    ctype = ["C"] * n + ["B"] * (nGpr+nAux+nGen)

    # options set
    OPTIONS = dict()
    '''
    OPTIONS['Display'] = 'off'
    OPTIONS['MaxTime'] = 100
    OPTIONS['IntegerTolerance'] = 10^(-9)
    OPTIONS['NumSol'] = numMultiSols
    OPTIONS['SearchMode'] = 1
    OPTIONS["limits/soluions"] = 10
    '''

    # time limit
    x_target = 0
    knockout = np.array([])
    tStart = time.time()
    for i in range(maxLoop): 
        # time limit
        timerun = time.time() - tStart
        if timerun > timeLimit:
            break

        # set ratio for each loop
        alpha = alpha + gap
        Aeqm[row_ratio, biomassID] = -1 * alpha
        Aeqm[row_ratio, targetID] = 1

        # solve milp
        optres, index = solver(fobj, A=Am, b=bm, Aeq=Aeqm, beq=beqm, lb=lb, ub=ub, vtype=ctype, opt="milp", pool=10)
        if optres.getStatus() != "optimal":
            continue
        pool = optres.getSols()

        # verify
        if not pool:
            continue
        knockout = np.zeros((nGen, 1))
        for poolSol in pool:
            x = [optres.getSolVal(poolSol, value) for value in index]
            if len(x) != n+nGpr+nAux+nGen:
                continue
            x_target, knockout_gene = verifyGeneKnock(model, x[n+nGpr+nAux:n+nGpr+nAux+nGen], biomassID, targetID, model.lb[biomassID, 0])
            if x_target > 0:
                knockout[:, 0] = knockout_gene
                return x_target, alpha, knockout

    return x_target, alpha, knockout

def maxParseGPR(model):
    # data scale
    nRxn = model.rxns.shape[0]
    nGen = model.genes.shape[0]
    nAux = 1
    parInfo = [0] * nRxn
    nRelation = 0
    nEqual = 0
    gprLabel = np.ones(nRxn)

    # get GPR rules
    for i in range(nRxn):
        gpr = model.grRules[i, 0] # string 
        if not gpr: # empty
            parInfoUnit = np.full((1, 3), "", dtype=object)
            gprLabel[i] = 0
            continue
        gpr = gpr[0]

        if ' ' in gpr: # multiple
            gprStr = "(" + gpr + ")"
            maxUnit = gprStr.count("(")
            nRelation = nRelation + 2 * maxUnit
            parInfoUnit = np.empty((maxUnit, 3), dtype=object)
            unitLabel = 1
            while unitLabel <= maxUnit:
                lbracket = [j for j in range(len(gprStr)) if gprStr[j] == "("]
                rbracket = [j for j in range(len(gprStr)) if gprStr[j] == ")"]
                rpoint = rbracket[0]
                for id in lbracket:
                    if id < rpoint:
                        lpoint = id
                    else:
                        break
                gprUnit = gprStr[lpoint+1:rpoint]
                gprStr = gprStr[0:lpoint] + "aux" + str(nAux) + gprStr[rpoint+1:]

                # assign one to one gene
                if unitLabel == maxUnit:
                    parInfoUnit[unitLabel - 1, 0] = "real" + str(i + 1)
                else:
                    parInfoUnit[unitLabel - 1, 0] = "aux" + str(nAux)
                    nAux += 1
                
                # single and/or boolean
                gprUnitElement = gprUnit.split(" ")
                parInfoUnit[unitLabel - 1, 1] = gprUnitElement[1]
                parInfoUnit[unitLabel - 1, 2] = gprUnitElement[0::2]
                unitLabel += 1
        else: # single
            parInfoUnit = np.empty((1, 3), dtype=object)
            parInfoUnit[0, 0] = "real" + str(i + 1)
            parInfoUnit[0, 1] = ""
            parInfoUnit[0, 2] = gpr
            nEqual += 1
        parInfo[i] = parInfoUnit
    nAux -= 1
    return parInfo, nRxn, nGen, nAux, nRelation, nEqual, gprLabel


def constructMatrix(model):
    # obtain parsing GPR info
    parInfo, nRxn, nGen, nAux, nRelation, nEqual, gprLabel = maxParseGPR(model)

    # get exist gene-rxn info 
    indGPR = [i for i in range(len(gprLabel)) if gprLabel[i] == 1]  # index of rxn has grRules
    indRea = [None] * len(gprLabel)
    nGpr = len(indGPR)  # number of rxn has grRules
    for i in range(nGpr):
        indRea[indGPR[i]] = i + 1 # rxn-index in list that trim rxns having no grRules

    # init matrix
    nMet = model.S.shape[0]
    # and/or complex relation
    gprMatrixLeft = np.zeros((nRelation,nRxn+nGpr+nGen+nAux))
    gprMatrixRight = np.zeros((nRelation,1))

    # lb*y <= v <= ub*y
    gprLabelDiag = np.diag(gprLabel)
    grCorrelationRxn = gprLabelDiag[indGPR, :]
    boundMatrixLeftTop = np.hstack((grCorrelationRxn,
                                      -1 * np.diag(model.ub[indGPR, :].flatten()),
                                      np.zeros((nGpr,nGen+nAux))))
    boundMatrixLeftBottom = np.hstack((-1 * grCorrelationRxn,
                                      np.diag(model.lb[indGPR, :].flatten()),
                                      np.zeros((nGpr,nGen+nAux))))
    boundMatrixLeft = np.concatenate((boundMatrixLeftTop, boundMatrixLeftBottom), axis=0)
    boundMatrixRight = np.zeros((2*nGpr,1))

    # one gene control
    gprEqualLeft = np.zeros((nEqual,nRxn+nGpr+nGen+nAux))
    gprEqualRight = np.zeros((nEqual,1))

    # S*v=0
    fluxEqualLeft = np.concatenate((model.S, np.zeros((nMet,nGpr+nGen+nAux))), axis=1)
    fluxEqualRight = np.zeros((nMet,1))

    # row point when add GPR relation
    rowPoint1 = 1
    rowPoint2 = 1

    # construct GPR matrix
    for i in range(len(parInfo)):        
        # get each reaction GPR parsing info 
        parInfoUnit = parInfo[i]  
        # has GPR
        if isinstance(parInfoUnit, int):
            continue
        maxUnit = parInfoUnit.shape[0]           
        # single gene control
        if maxUnit == 1 and not parInfoUnit[0, 1]:
            geneID = np.where(model.genes == parInfoUnit[0, 2])[0][0]
            gprEqualLeft[rowPoint1 - 1, nRxn + indRea[i] - 1] = 1
            gprEqualLeft[rowPoint1 - 1, nRxn + nGpr + nAux + geneID] = -1
            rowPoint1 += 1  
        # multiple gene control
        else:
            # iteratly get each auxilary GPR
            for j in range(maxUnit):
                eachAuxGPR = parInfoUnit[j, :]
                logicRelation = eachAuxGPR[1]
                n = len(eachAuxGPR[2])
                eachAuxGPR_ctrlGene = eachAuxGPR[2]
                if logicRelation == 'and':
                    # obj gene
                    objLoc = locateGene(nRxn, nGpr, nAux, indRea, eachAuxGPR[0], model)
                    gprMatrixLeft[rowPoint2 - 1, objLoc] = -1
                    gprMatrixLeft[rowPoint2, objLoc] = n
                    gprMatrixRight[rowPoint2 - 1, 0] = n - 1
                    # ctrl gene
                    for k in range(n):
                        objLoc = locateGene(nRxn, nGpr, nAux, indRea, eachAuxGPR_ctrlGene[k], model)
                        gprMatrixLeft[rowPoint2 - 1, objLoc] = 1
                        gprMatrixLeft[rowPoint2, objLoc] = -1
                    rowPoint2 += 2
                elif logicRelation == 'or':
                    # obj gene
                    objLoc = locateGene(nRxn, nGpr, nAux, indRea, eachAuxGPR[0], model)
                    gprMatrixLeft[rowPoint2 - 1, objLoc] = -n
                    gprMatrixLeft[rowPoint2, objLoc] = 1
                    # ctrl gene
                    for k in range(n):
                        objLoc = locateGene(nRxn, nGpr, nAux, indRea, eachAuxGPR_ctrlGene[k], model)
                        gprMatrixLeft[rowPoint2 - 1, objLoc] = 1
                        gprMatrixLeft[rowPoint2, objLoc] = -1
                    rowPoint2 += 2
                else:
                    raise ValueError('parsing GPR error')
                
    # construct matrix
    lessMatrixLeft = np.concatenate((boundMatrixLeft, gprMatrixLeft), axis=0)
    lessMatrixRight = np.concatenate((boundMatrixRight, gprMatrixRight), axis=0)
    equalMatrixLeft = np.concatenate((fluxEqualLeft, gprEqualLeft), axis=0)
    equalMatrixRight = np.concatenate((fluxEqualRight, gprEqualRight), axis=0)
    return lessMatrixLeft, lessMatrixRight, equalMatrixLeft, equalMatrixRight, nGpr, nAux, nGen, indGPR


def locateGene(nRxn, nGpr, nAux, indRea, objGene, model):
    geneNameLen = len(objGene)
    if geneNameLen > 2:
        geneNameHead = objGene[0:3]
        if geneNameHead == 'rea':
            ind = indRea[int(objGene[4:geneNameLen]) - 1]
            objLoc = ind + nRxn - 1 
        elif geneNameHead == 'aux':
            objLoc = int(objGene[3:geneNameLen]) + nRxn + nGpr - 1
        else:
            objLoc = np.where(model.genes == objGene)[0][0] + nRxn + nGpr + nAux   
    else:   
        objLoc = np.where(model.genes == objGene)[0][0] + nRxn + nGpr + nAux
    return objLoc

def verifyRatioGene(model, geneValue):
    # get gr rules and genes
    grRules_ori = model.grRules
    genes = model.genes
    num_gene = geneValue.shape[0]
    num_rxn = grRules_ori.shape[0]
    geneKnock = np.zeros((num_rxn, 1))

    # exclude array([], dtype='<U1')
    grRules = np.full((num_rxn, 1), "", dtype=object)
    for i in range(num_rxn):
        try:
            grRules[i][0] = grRules_ori[i][0][0]
        except IndexError:
            pass

    # replace and/or with */+
    grRules = grRules.astype(str)
    grRules = np.char.replace(grRules, 'or', '+')
    grRules = np.char.replace(grRules, 'and', '*')

    # replace gene to the value
    genesList = np.array([s[0][0] for s in genes])
    index = np.argsort(genesList)[::-1]
    genes1 = np.sort(genesList)[::-1] # some gene has similar name like YER060W_A and YER060W, YER060W_A need to be repleced first
    for i in range(num_gene):
        grRules = np.char.replace(grRules, genes1[i], str(geneValue[index[i], 0])) 
    for j in range(num_rxn):
        gpr = grRules[j, 0]
        if gpr == "":
            continue
        if eval(gpr) < 0.9:
            geneKnock[j, 0] = 1
    return geneKnock

def verifyGeneKnock(model, input_x, biomassID, targetID, min_bound):
    XX = 0
    length_input = len(input_x)
    input_xg = np.zeros((length_input, 1))
    for i in range(length_input):
        if input_x[i] > 0.1:
            input_xg[i, 0] = 1 
    var_x = verifyRatioGene(model,input_xg)
    index = var_x[:, 0] == 1
    model.lb[index, 0] = 0
    model.ub[index, 0] = 0
    model.lb[biomassID, 0] = 0

    optres, _ = solver(-model.c, Aeq=model.S, beq=model.b, lb=model.lb, ub=model.ub, opt="lp")
    EXITFLAG = optres.getStatus()
    if EXITFLAG != "optimal":
        return XX, input_xg
    FVAL = optres.getObjVal()
    if -1 * FVAL >= min_bound:
        # worst case
        model2 = copy.deepcopy(model)
        model2.lb[biomassID]=(-1)*FVAL
        model2.ub[biomassID]=(-1)*FVAL
        model2.c[biomassID]=0
        model2.c[targetID]=1
        optres2, index2 = solver(model2.c, Aeq=model2.S, beq=model2.b, lb=model2.lb, ub=model2.ub, opt="lp")
        x2 = optres2.getSols()[0]
        if x2[index2[targetID]] > 0.001:
            XX = x2[index2[targetID]]
    
    return XX, input_xg


def verifyRxn(model, knot, biomassID, targeID, LBbiomass):
    x_target = 0
    knot = np.abs(knot)
    delist = np.where(knot < 0.0000001)
    model.lb[delist[0], :] = 0
    model.ub[delist[0], :] = 0
    model.lb[biomassID][0] = 0

    optres, index = solver(-model.c, Aeq=model.S, beq=model.b, lb=model.lb, ub=model.ub, opt="lp")
    x = optres.getSols()[0]
    if x[index[biomassID]] > LBbiomass:
        x_target = x[index[targeID]]

    return x_target, delist[0]

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
    model = readm[modelName][0, 0]
    
    # get target ID
    target = args.targeMet
    try:
        targetName = np.where(model.mets == target)[0, 0]
    except IndexError:
        raise NameError("Check the target name, should be inside of model.mets")
    
    # get biomass ID
    biomassName = args.biomass
    if not biomassName:
        biomassID = np.nonzero(model.c)[0, 0]
    else:
        try:
            biomassID = np.where(model.rxns == biomassName)[0, 0]
        except IndexError:
            raise NameError("Check the biomass name, should be inside of model.rxns")
        
    # get carbon ID and oxygen ID
    carbon = args.carbon
    try:
        carbonID = np.where(model.rxns == carbon)[0, 0]
    except IndexError:
        raise NameError("Check the carbon name, should be an exchange reaction inside of model.rxns")
    oxygen = args.oxygen
    try:
        oxygenID = np.where(model.rxns == oxygen)[0, 0]
    except IndexError:
        raise NameError("Check the oxygen name, should be an exchange reaction inside of model.rxns")
    
    # assigen lbs to biomass carbon oxygen
    model.lb[biomassID, 0] = args.LBbiomass
    model.lb[carbonID, 0] = args.LBcarbon # <0
    model.lb[oxygenID, 0] = args.LBoxygen
    model.ub[carbonID, 0] = 0
    model.ub[oxygenID, 0] = 0

    # check target exchange rxns exist or not
    model, targetID, TMPR = introExchange(model,biomassID,[carbonID, oxygenID],targetName)
    if TMPR <= 0:
        raise ValueError("Theoritcal Maximum Production Rate is not positive.")
    # RatMethod for genes or rxns
    gap = TMPR/(args.maxLoop*args.LBbiomass)
    opt, _ = solver(-model.c, Aeq=model.S, beq=model.b, lb=model.lb, ub=model.ub, opt="lp")
    optSol = opt.getSols()[0]
    TMGR = -1 * opt.getSolObjVal(optSol)
    if args.type == "gene":
        x_obj, alpha, knockout = RatMethodGene(model, biomassID, targetID, TMGR, 10, gap, 1, 500)
    elif args.type == "reaction":
        x_obj, alpha, knockout = RatMethodRxn(model,biomassID,targetID,args.maxLoop,gap,args.timeLimit)
    else:
        raise NameError("Deletion type should be either \'gene\' or \'reaction\'")
    
    print(f"target production:{x_obj} \n alpha:{alpha}")