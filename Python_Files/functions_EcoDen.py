# ======================================================================================
# Copyright 2014  Swiss Federal Institute of Aquatic Science and Technology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#    
# Literature
# ==========
# Eggimann Sven, Truffer Bernhard, Maurer Max (2016): Economies of density for waste water. Water Research.
# 
                        
# Contact:   sven.eggimann@eawag.ch
# Version    1.0
# Date:      1.6.2015
# Autor:     Eggimann Sven
# ======================================================================================

import math, random, numpy, copy                    # Import Numpy for standard deviation
import numpy as np 

print("Sucseeful loading functions_transport")

def testVisualizeCatchement(scheduledRoute, InitialDepot, sinks):
    #----------------------
    #####TESET
    #----------------------
    testDistanceOneRound = 0
    distanceDepotCenter = 0
    totalTime = 0
 
    cnt = 0
    for i in scheduledRoute:
        if cnt == 1:
            distn = distanceCalc2dKmFactor(old, (i[1], i[2]), 1)       # Distance to clsest sink from current OST   
            totalTime += distn / float(vtruck) 
            testDistanceOneRound += distn
        old = (i[1], i[2])
        cnt = 1
                    
    # calc average distance to depot center
    for i in scheduledRoute:
        d = distanceCalc2dKmFactor(InitialDepot, (i[1], i[2]), 1)       # Distance to clsest sink from current OST   
        distanceDepotCenter += d
                

    print("-----")
    print("TEST CALCULATIONS Total Distance one round: " + str(testDistanceOneRound))
    print("Dist per OST: " + str(testDistanceOneRound/len(scheduledRoute)))
    print("TOTTIME: " + str(totalTime))
    print("Anzahl OST: " + str(len(scheduledRoute)))
    print("Aerage distance ode pot:  " + str(distanceDepotCenter/len(scheduledRoute)))
                
    import numpy as np
    import matplotlib.pyplot as plt
            
    xD = []; yD=[]
    for point in scheduledRoute:
        xD.append(point[1])
        yD.append(point[2])
    plt.scatter(xD,yD)
                
    print("InitialDepot:" + str(InitialDepot))
    depotX, depotY = InitialDepot[0], InitialDepot[1]  
    plt.scatter(depotX,depotY, s=80, color='red')      # Depot
                
    xSinks, ySinks = [],[]
        
    for i in sinks:
        xSinks.append(i[0])
        ySinks.append(i[1])
                
    plt.scatter(xSinks,ySinks, s=20, color='green')      # sinks
    plt.show() 

def roundTo10Steps(x, base=5):
    ''' Rounding Function to 10 steps
    
    Input:
    x       -    Number
    base    -    Rounding to 10 steps
    
    Output:
    rounded number
    
    '''
    return int(base * round(float(x) / base))
    
def getDistanceFactor(r, percentSA):
    '''
    depending on circle radius get distance factor
    
    Input:
    r                           -    Circle Radius
    percentSA                   -    Percentage for Sensitivity Analysis
    
    Output:
    fdScheduled                 -    Distance Factor scheduled
    fdUnscheduled               -    Distance Factor unscheduled
    
    '''
    fdScheduled = (0.0004 * r + 1.8536) * percentSA         # Derived distance weighting from case study
    fdUnscheduled = (-0.0005 * r + 1.2775) * percentSA      # Derived distance weighting from case study
    return fdScheduled, fdUnscheduled

def distanceCalc2dKmFactor(p0, p1, df):
    ''' Calculate 2d Distance'''

    distance = math.hypot(p0[0] - p1[0], p0[1] - p1[1]) 
    km = (float(distance) / 1000) * df 
    return km

def fromOldNodeToNewNode(currentOST, newOST, df, totTravelDistance, totTravelTime, hPersonellOperation, travelSpeed, notAlreadyAtnewNode):
    
    if notAlreadyAtnewNode == False:
        dist = distanceCalc2dKmFactor(currentOST, newOST, df)     # Distance from former position to new position
        travelTime = dist / float(travelSpeed)                                                  # Calculate Travel Time
        
        totTravelDistance += dist                                                               # Add Travel distance to total distance
        totTravelTime += travelTime
        hPersonellOperation += travelTime                                                      # Add Travel Time to current travel time

    return totTravelDistance, totTravelTime, hPersonellOperation 

def transformNetworkToDepot(streetVertices, depotCircle, streetNetworkCenter):
    '''
    transform street network in order that center is the circle center
    '''
    streetVerticesNew = []
            
    XShift = depotCircle[0] - streetNetworkCenter[0]
    YShift = depotCircle[1] - streetNetworkCenter[1]
            
    # Correct Vertices
    for i in streetVertices:
        newX = i[1] + XShift
        newY = i[2] + YShift
        z = [i[0], newX, newY, 0]
        streetVerticesNew.append(z)

    return streetVerticesNew

def checkIfAtHome(pemp):
    ''' check if anybody is at home'''
    randomDistr = random.uniform(0, 1)
    
    if randomDistr < pemp:
        atHome = True
    else:
        atHome = False
    
    return atHome
        
def emptyWWTPsInList(WWTPList, WWTP):
    '''
    This functions empties WWTPs in a list
    WWTPstoEmpty    -    List with WWTPs
    WWTP            -    WWTP to empty
    
    '''
    import copy
    WWTPEmpty = copy.copy(WWTP)
    
    # Emtpy WWTPs
    for i in WWTPList:
        for e in WWTPEmpty:
            if e[0] == i[0]:
                e[8] = 0
                break
    return WWTP

def controlA(WWTPstoEmpty):
    ''' control function'''
    #print(WWTPstoEmpty)
    for f in WWTPstoEmpty:
        if f[8] != 0:
            print("ERROR: NOT WWTP EMPTIED")
            prnt(".")               

def algorithmEvacuation(punkte, initialDepot, tpmax, sinks, travelSpeed, temptyTruck, df, cltruck, temptyOST, pemp):

    ''' This function calcululates the costs of emptying WWTPs with always checking if capacity reached and if yes travel to sink.  
    
    Input
    punkte             -    WWWTs in correct tour order
    initialDepot       -    Initial Depot
    tpmax              -    Max Working hours
    sinks              -    Deposits
    travelSpeed        -    Vehicle Travel Speed
    temptyTruck        -    Time needed to empty truck
    df   -    Distance Factor
    cltruck            -    Load Truck
    temptyOST          -    Time for emptying OST
    pemp               -    Emptying probability
    
    Output:
    totTravelDistance  -    Total travelled distance
    totTravelTime      -    Total travelled time
    totWorkingTime     -    Total working time
    
    '''   
    positionInWWTP, totTravelDistance, hPersonnel, totTravelTime, loadedOnTrucksOperation, totWorkingTime = 0, 0, 0, 0, 0, 0  # Parameters

    statistics_NrOfDEPOVISIT, statistics_Ausnahme = 0, 0
    
    visitedWWTP = []                                    # Visitied WWTPs
    currentOST = initialDepot                           # Initial copy
    emptyBecasueFarAwayOST = False
    while len(visitedWWTP) < len(punkte):               # As long not all WWTP have been visited
        #print("start " + str(emptyBecasueFarAwayOST))
        for i in punkte:
            if i[0] not in visitedWWTP:
                #print("NEW OST---------- " + str(i))
                #print("old OST---------  " + str(currentOST))
                #print("hPersonnel: " + str(hPersonnel))
                #print("lenVISIT: " + str(len(visitedWWTP)))
                newWWTPCoordinates = (i[1], i[2])       # Coordinates of WWTP to serve next
                notAlreadyAtnewNode = False   

                # Check if max work a day is reached. If yes, drive to depot and empty and return to new node the next day
                farAwayOST, hPersonnel, totTravelTime, totTravelDistance, loadedOnTrucksOperation, totWorkingTime, returnedToDepot, statistics_NrOfDEPOVISIT, statistics_Ausnahme = checkIfDayReached(statistics_Ausnahme, statistics_NrOfDEPOVISIT, emptyBecasueFarAwayOST, currentOST, newWWTPCoordinates, initialDepot, tpmax, loadedOnTrucksOperation, sinks, hPersonnel, totTravelTime, totTravelDistance, travelSpeed, totWorkingTime, temptyTruck, df, temptyOST)

                if returnedToDepot == True:             # 
                    #print("Return to depot..." + str(initialDepot))
                    currentOST = initialDepot           # Start iterating WWTPs again to test if any WWTP was missed because of somebody not beeing at home)
                    break                               # Start iterating path again
                
                # Check if OST emptying os possible. If no, truck still drives there
                atHome = checkIfAtHome(pemp)            # Check if anyone is at home
                if atHome == False:                     # Evacuation not possible
                    #print("not at home")
                    totTravelDistance, totTravelTime, hPersonnel = fromOldNodeToNewNode(currentOST, newWWTPCoordinates, df, totTravelDistance, totTravelTime, hPersonnel, travelSpeed, notAlreadyAtnewNode)
                    currentOST = (newWWTPCoordinates)   # Truck is not at OSt which cannot be emptied
                    continue
                #print(len(visitedWWTP))
                #print("hPersonel:                  " + str(hPersonnel))
                #print("totTravelTime:              " + str(totTravelTime))
                #print("totTravelDistance:          " + str(totTravelDistance))
                #print("totWorkingTime:             " + str(totWorkingTime))
                #print("loadedOnTrucksOperation:    " + str(loadedOnTrucksOperation))
                #print("returnedToDepot:            " + str(returnedToDepot))
                
                # Until current WWsTP is fully emptied    
                totTravelDistance, totTravelTime, hPersonnel = fromOldNodeToNewNode(currentOST, newWWTPCoordinates, df, totTravelDistance, totTravelTime, hPersonnel, travelSpeed, notAlreadyAtnewNode)
                    
                capacityReached, stillEmptyCapacity = checkIfthereIsSpaceForAnotherWWTP(loadedOnTrucksOperation, cltruck, i)            # Check if truck capacity is reached
                punkte, loadedOnTrucksOperation, hPersonnel, totWorkingTime = emptyWWTP(punkte, newWWTPCoordinates, stillEmptyCapacity, loadedOnTrucksOperation, temptyOST, totWorkingTime, hPersonnel)                    # Empty as much as possible in next WWTP. 
                #print(".." + str(totTravelDistance))
                #print(".." + str(stillEmptyCapacity))
                #print(".." + str(hPersonnel))
                if capacityReached == 1 or farAwayOST == True:  # Either WWTP not empty or farAwayOST          
                    #print("Truck full or farAwayOST: " + str(farAwayOST))
                    SinkCoordainte = getClosestSink(newWWTPCoordinates, sinks, df)                                         # Search closest sink
                    distanceToSink = distanceCalc2dKmFactor(newWWTPCoordinates, SinkCoordainte, df)                        # distance to Sink               
                    travelTime = distanceToSink / float(travelSpeed)                                                                     # Travel time to Sink
                        
                    totTravelTime += travelTime                                                                                          # Summen time
                    hPersonnel += travelTime                                                                                             # Summen current daily time
                    totTravelDistance += distanceToSink                                                                                  # summen total travelled distance

                    loadedOnTrucksOperation, totWorkingTime = emptyAtSink(loadedOnTrucksOperation, totWorkingTime, temptyTruck)          # Empty Sink    
                    currentIsEmpty = checkIfCurrentWWTPIsEmpty(punkte, i)                                                                # If new is empty, drive to new otherwise keep continung emptying same wwtp
                    #print("currentIsEmpty: " + str(currentIsEmpty))
                        

                    if farAwayOST == False:                     # If regular OST
                        currentOST = (SinkCoordainte)       # Current position
                        if currentIsEmpty == True:              # If regular OST is empty
                            #print("A")
                            _ = 0
                        else:
                            #print("B")
                            # Start iterating again (then a half-empty OST is first found
                            break           # start iteration again
                    else:

                        if currentIsEmpty == True: # WWTP is empty, return immediatly to Depot because it is a FarAway OST
                            #print("C")
                            # return immediately to Depot    
                            distanceToDepotOrig = distanceCalc2dKmFactor(SinkCoordainte, initialDepot, df)                        # distance to Sink               
                            travelTimeToDepotOrig = distanceToDepotOrig / float(travelSpeed)                                                                     # Travel time to Sink
                                
                            totTravelTime += travelTimeToDepotOrig                                                                                          # Summen time
                            hPersonnel += travelTimeToDepotOrig                                                                                             # Summen current daily time
                            totTravelDistance += distanceToDepotOrig                                                                                  # summen total travelled distance
                            currentOST = (initialDepot)   # Current position
                            #print("totTravelTime:      " + str(totTravelTime))
                            #print("hPersonnel:         " + str(hPersonnel))
                            #print("totTravelDistance:  " + str(totTravelDistance))
                            #print("currentOST:         " + str(currentOST))
                                
                            hPersonnel = 0 # Set to zero again
                            emptyBecasueFarAwayOST = False
                        else: # Continue emptying until all is empties
                            #print("D")
                            emptyBecasueFarAwayOST = True
                            currentOST = (SinkCoordainte)   # Current position
                            break # start iterating again

                else:                                   # Capacity not reach in truck
                    currentIsEmpty = True               # If truck is not full after the whole tank was emptied, the OSt is certainly empty. Then move to next OST
                    currentOST = (newWWTPCoordinates)   # punkte last visited    

                # Return to depot a the very end of iteration
                if positionInWWTP == len(punkte)-1:
                        #print("Final return to the depot")
         
                        # Search closest sink, empty tank
                        SinkCoordainte = getClosestSink(newWWTPCoordinates, sinks, df)               # Drive from new to sink
                        dist = distanceCalc2dKmFactor(newWWTPCoordinates, SinkCoordainte, df)                
                        travelTime = dist / float(travelSpeed)
                        
                        totTravelTime += travelTime
                        hPersonnel += travelTime
                        totTravelDistance += dist
                        
                        loadedOnTrucksOperation, totWorkingTime = emptyAtSink(loadedOnTrucksOperation, totWorkingTime, temptyTruck)    # Empty Sink 
                        
                        # Drive from Sink to depot
                        dist = distanceCalc2dKmFactor(SinkCoordainte, initialDepot, df)
                        totTravelDistance += dist
                        travelTime = dist / float(travelSpeed)
                        totTravelTime += travelTime
                    
                #print("Add to visited")
                positionInWWTP += 1
                visitedWWTP.append(i[0])
                    
                #print("ZUWACHS: " + str(totTravelDistance - totTravelDistanceOld))
    #print("===========================================================")
    #print("STATS ROUTE")
    #print("statistics_Ausnahme:         " + str(statistics_Ausnahme))
    #print("statistics_NrOfDEPOVISIT: " + str(statistics_NrOfDEPOVISIT))
    #print("LEN PUNKTdddE: " + str(punkte))

    return totTravelDistance, totTravelTime, totWorkingTime

def algorithmServiceTour(punkte, Depot, InitialDepot, tpmax, travelSpeed, df, tservice):

    ''' This function calculates costs of travelling to a OST an performing tasks
    
    punkte                    -    OST in correct order
    Depot                     -    Depot
    InitialDepot              -    InitialDepot
    tpmax                     -    maxium working time
    sinks                
    travelSpeed
    df
    tservice                  -    Service time at WWTP
    
    Output:
    totTravelDistance         -    Tot travel distance
    totTravelTime             -    Tot travel time
    totWorkingTime            -    tot working time
    '''
    positionInWWTP, totTravelDistance, totTravelTime, totWorkingTime, hPersonell = 0, 0, 0, 0, 0        # Parameters

    for i in punkte:
        newWWTPCoordinates = (i[1], i[2])                                                                # Coordinates of WWTP to serve next

        # Check if max work a day is reached. If yes, drive to depot and empty and return to new node the next day
        hPersonell, totWorkingTime, totTravelTime, totTravelDistance, returnToDepot = checkServiceTimeConstraint(Depot, newWWTPCoordinates, InitialDepot, tpmax, hPersonell, totTravelTime, totTravelDistance, travelSpeed, df, totWorkingTime, tservice)

        if returnToDepot == True:
            Depot = InitialDepot
        else:
            Depot = (newWWTPCoordinates)  # punkte last visited    

        # Return to depot a the very end of iteration
        if positionInWWTP == len(punkte) - 1:
            dist = distanceCalc2dKmFactor(newWWTPCoordinates, InitialDepot, df)
            totTravelDistance += dist
            travelTime = dist / float(travelSpeed)
            totTravelTime += travelTime
            
        positionInWWTP += 1
    return totTravelDistance, totTravelTime, totWorkingTime

def clarkeAndWright(WWTPInput, InitialdepotOperation, df):
    '''
    
    Input:
    
    Output:
    
    
    '''
    #print("start creating clakre & wright" + str(len(WWTPInput)))
    #print("----")
    WWTP = copy.deepcopy(WWTPInput)
    
    ID_INITALDEPOT, ID_INITALDEPOT_orig = -1, -1
    depotClarke = (ID_INITALDEPOT, InitialdepotOperation[0], InitialdepotOperation[1])  # Clark Depot
    nrOfConnectionNodes, path, connectedNodes = len(WWTP), {}, []

    distanceMatrix = calculateDistanceMatrix(WWTP, depotClarke, df)       # Calculate distance matrix 
    #print("Distance Matrix is calculated: ")
    #print("DEBOP: " + str(depotClarke))

    savingList = calculateSavings(distanceMatrix, depotClarke, WWTP)                    # Calculate savings
    
    #print("FF: " + str(savingList[:10]))
    #print("Savingas are calculated: " + str(len(savingList)))
 
    savingList, FROMNODE, TONODE = getMostSaving(savingList)                            # Connect initial triangle. Find node with most savings

    dFROMNODE = distanceMatrix[depotClarke[0]][FROMNODE]
    dTONODE = distanceMatrix[depotClarke[0]][TONODE]
    DINBETWEEN = distanceMatrix[FROMNODE][TONODE]
        
    path[depotClarke[0]] = {FROMNODE: dFROMNODE}
    path[FROMNODE] = {TONODE: DINBETWEEN}
    path[TONODE] = {depotClarke[0]: dTONODE}
        
    # new and hopefully faster
    savingList.sort()
    savingList = savingList[::-1]

    # Connect the rest of the points
    if len(WWTPInput) > 2: # Otherwise only one option (we do not consider direction)
        path = getMostSavingInPathDIRECT(nrOfConnectionNodes, savingList, path, distanceMatrix, depotClarke, connectedNodes)  # Generates saving edge

    WWTsavingListEW = []

    new = 0
    while len(WWTsavingListEW) != len(WWTPInput): #new != ID_INITALDEPOT_orig:
        toID = path[ID_INITALDEPOT]
        
        for e in toID:
            new = e
            break
        if new == ID_INITALDEPOT_orig:
            # Get wwtp
            for e in WWTP:
                if e[0] == ID_INITALDEPOT and ID_INITALDEPOT != ID_INITALDEPOT_orig:
                    WWTsavingListEW.append(e)
                    break
            break
        else:
            # Get wwtp
            for e in WWTP:
                if e[0] == new and new != ID_INITALDEPOT_orig:
                    WWTsavingListEW.append(e)
                    break
            ID_INITALDEPOT = new
    return WWTsavingListEW                


def calculateDistanceMatrix(WWTP, depot, df):
    '''Calculate distance matrix. 2D straight'''
    WWTP.insert(0, depot)  # Add depot
    distanceMatrix, cnt = {}, 0  # Distance matrix     

    for fromNode in WWTP:
        cnt += 1                        
        distanceMatrix[fromNode[0]] = {}  # Add empty dictionary
        
        # En Reched of WWTP (no distance aclualtionp ossibel
        if fromNode[0] == WWTP[-1][0]:
            return distanceMatrix
        
        for toNode in WWTP[cnt:]:  # Calculate distances with djikstra between all pairs
            distance = distanceCalc2dKmFactor((fromNode[1], fromNode[2]), (toNode[1], toNode[2]), df)
            distanceMatrix[fromNode[0]][toNode[0]] = distance


def calculateSavings(distanceMatrix, depot, nodeList):
    savingList, cnt = [], 1

    for fromNode in nodeList:
        if fromNode[0] == nodeList[-1][0]:  # If last elements
            return savingList
        for ToNode in nodeList[cnt:]:
            if fromNode[0] != ToNode[0] and fromNode[0] != depot[0]:

                d_depotToSecondNode = distanceMatrix[depot[0]][fromNode[0]]
                d_depotToFirstNode = distanceMatrix[depot[0]][ToNode[0]]
                distancebetween = distanceMatrix[fromNode[0]][ToNode[0]]
                savings = d_depotToFirstNode + d_depotToSecondNode - distancebetween
                #print("---")
                #print(d_depotToSecondNode)
                #print(d_depotToFirstNode)
                #print(distancebetween)
                #print(savings)
                
                savingList.append((savings, fromNode[0], ToNode[0], d_depotToFirstNode, d_depotToSecondNode, distancebetween))
        cnt += 1
    return savingList

def getMostSaving(savingList):
    minDit = 0  # scrapdistance
    zahler = -1
    for i in savingList:
        zahler += 1
        if i[0] > minDit:      
            minDit, FROMNODE, TONODE = i[0], i[1], i[2]  # shortest distance      # the flow flows to this node.  # the flow starts here and flows to fromnode      
            deletPosition = zahler
            #print("M: " + str(minDit))
    del savingList[deletPosition]

    return savingList, FROMNODE, TONODE

def getMostSavingInPathDIRECT(nrOfConnectionNodes, savingList, P, distanceMatrix, depot, connectedIntermediatenodes):
    
    initialCopysavingList = list(savingList)
    depotCoordinate = depot[0]
    #while nrOfConnectionNodes >= len(P):
    
    #while nrOfConnectionNodes + 1 > len(P): # Added one because of depot
    #while initialCopysavingList > 0: 
    #print("A: " + str(nrOfConnectionNodes))  
    
    while nrOfConnectionNodes + 1 != len(P):
        #print("FORTSCHRITT: " + str(len(P)))
        #print("len initi:   " + str(len(initialCopysavingList)))
        #print(len(P))
        #print("..")
        #print("SAVING LIST: " + str(savingList))
        #if len(savingList) == 0:
        if len(initialCopysavingList) == 0 or len(savingList) == 0:    
            print("ERROR with saving list: ")
            prnt("........")
        #print(nrOfConnectionNodes)
        #print("===========================")
        
        # Remove interconnecting connection from savingList
        savingListNew = []
        for i in initialCopysavingList:
            if i[1] in P and i[2] in P:
                _ = 0
            else:
                savingListNew.append(i)
        #initialCopysavingList = copy.deepcopy(savingListNew)
        initialCopysavingList = list(savingListNew)

        #savingList = copy.deepcopy(initialCopysavingList)  
        savingList = list(initialCopysavingList) 
        #print("AFTER I: " + str(savingList))  
  
        # Iterate a long all pairs have been checked or a new connection took place
        ext = 1
        deletPositionSavingList = -1
        while len(savingList) > 0 and ext == 1:
            deletPositionSavingList += 1
            
            #print("savingList[0]: " + str(savingList[0]))
            FROMNODE, TONODE = savingList[0][1], savingList[0][2]  # Get shortest distance
            deletPosition = 0  

            del savingList[deletPosition]  # Remove connection
            
            # get possible nodes where a connection is possible
            for i in P[depotCoordinate]:
                fromDepot = i
                break
            
            for i in P:
                for e in P[i]:
                    if e == depotCoordinate:
                        toDepot = i
                        break
            #print("Endpunkte:   " + str(toDepot) + "   " + str(fromDepot))
            #print("FROMNODE: " + str(FROMNODE) + "   " + str(TONODE))
            #print("connectedIntermediatenodes: " + str(connectedIntermediatenodes))
            #print("----------------------------------------------------------------")
            
            # If connetion between WWTP already exists
            if FROMNODE == toDepot and TONODE == fromDepot or FROMNODE == fromDepot and TONODE == toDepot:  # If not connecting possible nodes which can be connected
                del initialCopysavingList[deletPosition]  # Remove connection
            else:
                # Test if end nodes are next to new shortest saving path nodes
                if FROMNODE == fromDepot or TONODE == fromDepot or FROMNODE == toDepot or TONODE == toDepot: 
                    if FROMNODE == fromDepot or TONODE == fromDepot:
                        if FROMNODE == fromDepot:
                            if TONODE in connectedIntermediatenodes:            # and FROMNODE in connectedIntermediatenodes:
                                continue                                        # Return to while                              # Select next coordinate pair
                            dist = distanceMatrix[depotCoordinate][TONODE]      # find distances
                            dist2 = distanceMatrix[FROMNODE][TONODE]            # find distances        
                            P[depotCoordinate] = {TONODE: dist}                 # Change P
                            P[TONODE] = {FROMNODE: dist2}                       # Change P
                            connectedIntermediatenodes.append(fromDepot)        # List with already connected nodes
                        else:
                            if FROMNODE in connectedIntermediatenodes:          # and FROMNODE in connectedIntermediatenodes:
                                continue                                        # Select next coordinate pair
                            dist = distanceMatrix[depotCoordinate][FROMNODE]    # find distances
                            dist2 = distanceMatrix[FROMNODE][fromDepot]         # find distances    
                            P[depotCoordinate] = {FROMNODE: dist}               # Change P
                            P[FROMNODE] = {fromDepot: dist2}                    # Change P
                            connectedIntermediatenodes.append(fromDepot)        # List with already connected nodes
                    if FROMNODE == toDepot or TONODE == toDepot:    
                        if toDepot == FROMNODE:
                            if TONODE in connectedIntermediatenodes:            # and FROMNODE in connectedIntermediatenodes:
                                continue                                        # Select next coordinate pair
                            dist = distanceMatrix[toDepot][TONODE]              # find distances
                            dist2 = distanceMatrix[depotCoordinate][TONODE]            # find distances
                            P[toDepot] = {TONODE: dist}                         # Change P
                            P[TONODE] = {depotCoordinate: dist2}                # Change P
                            connectedIntermediatenodes.append(toDepot)          # List with already connected nodes
                        else:
                            if FROMNODE in connectedIntermediatenodes:          # and FROMNODE in connectedIntermediatenodes:
                                continue
                            dist = distanceMatrix[FROMNODE][toDepot]            # find distances
                            dist2 = distanceMatrix[depotCoordinate][toDepot]    # find distances
                            P[toDepot] = {FROMNODE: dist}                       # Change P
                            P[FROMNODE] = {depotCoordinate: dist2}              # Change P
                            connectedIntermediatenodes.append(toDepot)          # List with already connected nodes
                    del initialCopysavingList[deletPositionSavingList]          # Remove connection
                    ext = 0
    return P

def distanceCalc2dKm(p0, p1):
    ''' Calculate 2d Distance'''

    distance = math.hypot(p0[0] - p1[0], p0[1] - p1[1]) 
    km = (float(distance) / 1000) 
    return km

def averageNearestNeighborClustering(buildings, areaCircle):
    '''
    input: 
    buildings    -    Buildings
    areaCircle     -    Area of circle
    
    Output:
    ANN          -    ANN    
    '''
    A = areaCircle * 1000000  # [m2]                                                   
    nrNodes = len(buildings)            

    
    # Calculate nearest feature for each point
    SumOfNearestNeigbhour = 0
    for b in buildings:
        
        # Window Search
        nearDi = 99999999999999  # Nearest Distance of a point
        for e in buildings:  # Search distance to closest point
            dst = distanceCalc2dKm((b[1], b[2]), (e[1], e[2]))
            if dst != 0 and dst < nearDi:
                nearDi = dst

        SumOfNearestNeigbhour += nearDi

    SumOfNearestNeigbhour = SumOfNearestNeigbhour * 1000  # [m] Needed in Meters
    
    De = 0.5 / ((nrNodes / float(A)) ** 0.5)  # 
    Do = SumOfNearestNeigbhour / nrNodes  # distance in m                
    ANN = float(Do) / float(De)  # ANN
    # print("===========================")
    # print("buildings:      " + str(len(buildings)))
    # print("De:             " + str(De))
    # print("Do:             " + str(Do))
    # print("A:                 " + str(areaCircle))
    # print("CLUSTERING:     " + str(ANN))
    return ANN

def checkIfIsToocloseToBuilding(removedPntList, minInterWWTPDistance, random_X, random_Y):
    for f in removedPntList:
        innerWindowX_max = f[0] + minInterWWTPDistance * 1000 # [m]
        innerWindowX_min = f[0] - minInterWWTPDistance * 1000 # [m]
        innerWindowY_max = f[1] + minInterWWTPDistance * 1000 # [m] 
        innerWindowY_min = f[1] - minInterWWTPDistance * 1000 # [m]
        
        if random_X < innerWindowX_max and random_X > innerWindowX_min and random_Y < innerWindowY_max and random_Y > innerWindowY_min:  # KORRIGENDA OR ZU UND
            isCloseToBuilding = True
            break
        else:
            isCloseToBuilding = False
    return isCloseToBuilding

def generateSinks(extent, minInterSinkDistance, r, SinkDensity):
    '''
    # This function creates random sinks depending on the density
    
    Input: 
    extent                            -    Extent [[topLeftPoint], [topRightPoint]]
    minInterSinkDistance              -    Minimum Distance between sinks
    r                      -    Radius of Circle
    SinkDensity                       -    Nr of Sinks per km2
    
    return
    sinkList                              -    Sinks
    circleAre                         -    Area of Cirlse [km2]
    '''
    sinkList = []
    TopLefXGlobal, TopLefYGlobal, BottomRightXGlobal, BottomRightYGlobal = extent[0][0], extent[1][1], extent[1][0], extent[0][1]  # Get Coordinates of Extent
    
    # Center of Circle
    centerExtentX = TopLefXGlobal + ((BottomRightXGlobal - TopLefXGlobal) / 2)
    centerExtentY = BottomRightYGlobal + ((TopLefYGlobal - BottomRightYGlobal) / 2)
    areaCircle = math.pi * pow(r, 2)                                     # Circle Area
    nrOfOST = int(SinkDensity * areaCircle)                                        # Number of Points

    # Add Distance in order that circle
    TopLefXGlobal = centerExtentX - r * 1000                             
    TopLefYGlobal = centerExtentY + r * 1000
    BottomRightXGlobal = centerExtentX + r * 1000 
    BottomRightYGlobal = centerExtentY - r * 1000 
 
    # Create Random sinkList
    while len(sinkList) < nrOfOST:
        random_X = random.uniform(TopLefXGlobal, BottomRightXGlobal)  # Select random X
        random_Y = random.uniform(BottomRightYGlobal, TopLefYGlobal)  # Selext random Y
        randomPnt = [random_X, random_Y]                              # Random Point
        
        distanceFromMiddleOfExtent = distanceCalc2dKm((centerExtentX, centerExtentY), (random_X, random_Y)) 
        
        if distanceFromMiddleOfExtent < r:
            
            # Check if minimum distance criteria is fulfilled
            if len(sinkList) > 0:
                dCrit = True
                for i in sinkList:
                    distance = distanceCalc2dKm([random_X, random_Y], [i[0], i[1]])

                    if distance <= minInterSinkDistance:  
                        dCrit = False
                if dCrit == True:
                    sinkList.append(randomPnt)
            else:
                sinkList.append(randomPnt)
    
    # If not a single sink was created, create sink at center
    if len(sinkList) == 0:
        sinkList = [[centerExtentX, centerExtentY]]
    
    return sinkList, areaCircle
    
def checkIfDayReached(statistics_Ausnahme, statistics_NrOfDEPOVISIT, emptyBecasueFarAwayOST, currentOST, newWWTPCoordinates, Initialdepot, tpmax, loadedOnTrucks, sinks, hPersonell, totTravelTime, totTravelDistance, travelSpeed, totWorkingTime, temptyTruck, df, temptyOST): 
    ''' THis function checks if max work a day is reached. IF yes, add distance to depot and return to new node'''

    # Time for Depot to OST. Used to check if an exception is allowed
    initialDepotToNewOST = distanceCalc2dKmFactor(Initialdepot, newWWTPCoordinates, df)       # Distance to clsest sink from current OST            
    timeInitialDepotToNewOST = initialDepotToNewOST / float(travelSpeed)                                      # Travel 
        
    # Current OST to closest sink
    SinkCoordainte = getClosestSink(currentOST, sinks, df)                    # Closest sink from current OST
    distToSink = distanceCalc2dKmFactor(currentOST, SinkCoordainte, df)       # Distance to clsest sink from current OST            
    travelTimeCurrentOstToSINK = distToSink / float(travelSpeed)                                      # Travel tim from current OST to closest sink
    
    ## All Distances needed for new OST
    # from current Ost to new OST
    distNewToOST = distanceCalc2dKmFactor(currentOST, newWWTPCoordinates, df)                    # Distance from current WWTP to depot
    timeToNewOST = distNewToOST / float(travelSpeed)                          # Time neede from current OST to new OST

    # from new Ost to closest Deposit
    SinkCoordainteNewOST = getClosestSink(newWWTPCoordinates, sinks, df)                      # Closest sink from new OST
    distToSinkNewOST = distanceCalc2dKmFactor(newWWTPCoordinates, SinkCoordainteNewOST, df)   # Distance to clsest sink from new OST            
    travelTimeNewOST = distToSinkNewOST / travelSpeed                                                       # Travel tim from current OST to new sink
    
    # from closest sink of new OST to Depot
    distSinkToDepot = distanceCalc2dKmFactor(SinkCoordainteNewOST, Initialdepot, df)                    # Distance from current WWTP to depot
    travelTimeNewSinkToDepot = distSinkToDepot / float(travelSpeed)
    
    # Time from current OST
    totTimeNeeded = timeToNewOST + temptyOST + travelTimeNewOST + travelTimeNewSinkToDepot + temptyTruck # from current ot new OST + emptying time + time to closest sink + time from cloest sink to depot
    
    # Theoretical used to test for exception
    totTimeNeededTheoretical = timeInitialDepotToNewOST + temptyOST + travelTimeNewOST + travelTimeNewSinkToDepot + temptyTruck # from current ot new OST + emptying time + time to closest sink + time from cloest sink to depot
    
    #print("THEORETICAL TIME: " + str(totTimeNeededTheoretical))
    #print(statistics_Ausnahme)
    
    ### If overtime from visist of exceptionally far away node based on current situation is smaller than
    ### overtime in case we woulc start at a depot, vitis the exceptionally far away node.
    ### Otherwise return to depot.  
    
    ##print("theoreticalOverTimeIfStartedFromDepot: " + str(totTimeNeededTheoretical))
    ##print("OverTimeIfStartedFromCurrentWWTP:      " + str(totTimeNeeded))
    ##print(Initialdepot)
    ##print(newWWTPCoordinates)
    
    # In case a single OST is too far away to serve one day, make an exception
    if totTimeNeededTheoretical > tpmax or emptyBecasueFarAwayOST == True:
        statistics_Ausnahme += 1
        #print("Emptying far away plant from current situation which results in heavy overtime "  + str(tpmax))
        #print(totTimeNeededTheoretical)
        #print("travelSpeed: " + str(travelSpeed))
        #print("timeInitialDepotToNewOST      " + str(timeInitialDepotToNewOST))
        #print("temptyOST                     " + str(temptyOST))
        #print("travelTimeNewOST              " + str(travelTimeNewOST))
        #print("travelTimeNewSinkToDepot      " + str(travelTimeNewSinkToDepot))
        #print("temptyTruck                   " + str(temptyTruck))
        
        returnedToDepot = False # Only return to depot after farAwaysOST is emptied
        farAwayOST = True
        return farAwayOST, hPersonell, totTravelTime, totTravelDistance, loadedOnTrucks, totWorkingTime, returnedToDepot, statistics_NrOfDEPOVISIT, statistics_Ausnahme
        # Is more time efficient to visit from current situtation than from depot

    # Check if emptying is possible! As the w WWTP could not be reached make a tour via depot
    if totTimeNeeded + hPersonell >= tpmax:    # too much time is used
        #print("Emptying not possible" + str(totTimeNeeded + hPersonell))
        
        # Search closest sink
        totTravelTime += travelTimeCurrentOstToSINK
        totTravelDistance += distToSink

        # Set load on trucks to zero and add emptying truck to working time
        loadedOnTrucks, totWorkingTime = emptyAtSink(loadedOnTrucks, totWorkingTime, temptyTruck)    # Empty Sink 
        
        # Travel from Sink to Starting point
        distToDepot = distanceCalc2dKmFactor(SinkCoordainte, Initialdepot, df)                    # Distance from current WWTP to depot
        travelTime = distToDepot / float(travelSpeed)

        totTravelTime += travelTime
        totTravelDistance += distToDepot
        
        #print("Abweichung von h max: " + str(tpmax-(hPersonell+ travelTimeCurrentOstToSINK + temptyTruck + travelTime + temptyTruck)))
        hPersonell = 0                                                                                          # New day is started for personel 
    
        returnedToDepot = True                                                                                  # Return to Depot
        statistics_NrOfDEPOVISIT += 1
        farAwayOST = False
        return farAwayOST, hPersonell, totTravelTime, totTravelDistance, loadedOnTrucks, totWorkingTime, returnedToDepot, statistics_NrOfDEPOVISIT, statistics_Ausnahme#, ausnahme
    else:
        #print("Emptying  possible")
        returnedToDepot = False
        farAwayOST = False
        return farAwayOST, hPersonell, totTravelTime, totTravelDistance, loadedOnTrucks, totWorkingTime, returnedToDepot, statistics_NrOfDEPOVISIT, statistics_Ausnahme#, ausnahme

def checkServiceTimeConstraint(currentOST, newWWTPCoordinates, Initialdepot, tpmax, hPersonell, totTravelTime, totTravelDistance, travelSpeed, df, totWorkingTime, tservice):
    ''' THis function checks if max work a day is reached. IF yes, add distance to depot and return to new node'''
    
    # From current OST to next OST
    distToNewost = distanceCalc2dKmFactor(currentOST, newWWTPCoordinates, df)         # Distance from former position to new position
    travelTimeTOnewOST = distToNewost / float(travelSpeed)                                                      # Travel Time

    ## All Distances needed for new OST
    # From current Ost to new OST
    distNewToOST = distanceCalc2dKmFactor(currentOST, newWWTPCoordinates, df)                    # Distance from current WWTP to depot
    timeToNewOST = distNewToOST / float(travelSpeed)                          # Time neede from current OST to new OST

    # from next OST to Depot
    distNewOstToDepot = distanceCalc2dKmFactor(newWWTPCoordinates, Initialdepot, df)                    # Distance from current WWTP to depot
    travelTimeNewOstToDepot = distNewOstToDepot / float(travelSpeed)
    
    # Total potential needed time
    totTimeNeeded = timeToNewOST + tservice + travelTimeNewOstToDepot   # Time to drive to new OST + time spent at OST + time spent for driving to Depot

    # In case a single OST is too far away to serve one day, make an exception
    if totTimeNeeded > tpmax:
        returnToDepot = False
        return hPersonell, totWorkingTime, totTravelTime, totTravelDistance, returnToDepot
    
    # If max working time would be reached with new OST, return to Depot 
    if hPersonell + totTimeNeeded >= tpmax:  
        distToDepot = distanceCalc2dKmFactor(currentOST, Initialdepot, df)                    # Distance from current WWTP to depot
        dist = distToDepot #+ distFromDepotToNewNode

        travelTime = dist / float(travelSpeed)
        
        totTravelDistance += dist
        totTravelTime += travelTime
        
        hPersonell = 0  # New day is started for personel

        returnToDepot = True
        return hPersonell, totWorkingTime, totTravelTime, totTravelDistance, returnToDepot
    else:
        
        hPersonell += timeToNewOST + tservice
        
        # Go to next WWTP and perform service
        totTravelDistance += distToNewost
        totTravelTime += travelTimeTOnewOST
        totWorkingTime += tservice
        returnToDepot = False
        
        return hPersonell, totWorkingTime, totTravelTime, totTravelDistance, returnToDepot

def getCircleExtent(extent, r):
    '''
    Create random distribution within a circle within the extent for a given number of OST.
    
    Input: 
    mindistBetweenWWTPs               -    minimum distance between OST
    extent                            -    [[topLeftPoint], [topRightPoint]]
    nrOfOST                           -    Nr of OST
    r                                 -    radius
    '''

    randomPnts = []
    TopLefXGlobal, TopLefYGlobal, BottomRightXGlobal, BottomRightYGlobal = extent[0][0], extent[1][1], extent[1][0], extent[0][1]  # Get Coordinates of Extent

    centerExtentX = TopLefXGlobal + ((BottomRightXGlobal - TopLefXGlobal) / 2)
    centerExtentY = BottomRightYGlobal + ((TopLefYGlobal - BottomRightYGlobal) / 2)

    # Add Distance in order that circle
    TopLefXGlobal = centerExtentX - (r * 1000)         
    TopLefYGlobal = centerExtentY + (r * 1000)         
    BottomRightXGlobal = centerExtentX + (r * 1000)   
    BottomRightYGlobal = centerExtentY - (r * 1000)    
    
    circleExtent = [(TopLefXGlobal, BottomRightYGlobal), (BottomRightXGlobal, TopLefYGlobal)]
    return circleExtent

def createRandomDisributionCircle(mindistBetweenWWTPs, extent, nrOfOST, r):
    '''
    Create random distribution within a circle within the extent for a given number of OST.
    
    Input: 
    mindistBetweenWWTPs               -    minimum distance between OST
    extent                            -    [[topLeftPoint], [topRightPoint]]
    nrOfOST                           -    Nr of OST
    r                                 -    radius
    '''
    randomPnts = []
    TopLefXGlobal, TopLefYGlobal, BottomRightXGlobal, BottomRightYGlobal = extent[0][0], extent[1][1], extent[1][0], extent[0][1]  # Get Coordinates of Extent

    centerExtentX = TopLefXGlobal + ((BottomRightXGlobal - TopLefXGlobal) / 2)
    centerExtentY = BottomRightYGlobal + ((TopLefYGlobal - BottomRightYGlobal) / 2)

    # Add Distance in order that circle
    TopLefXGlobal = centerExtentX - (r * 1000)         
    TopLefYGlobal = centerExtentY + (r * 1000)         
    BottomRightXGlobal = centerExtentX + (r * 1000)   
    BottomRightYGlobal = centerExtentY - (r * 1000)    
    
    circleExtent = [(TopLefXGlobal, BottomRightYGlobal), (BottomRightXGlobal, TopLefYGlobal)]
    
    # Create Random randomPnts
    randomID = 0
    while len(randomPnts) < nrOfOST:
        random_X = random.uniform(TopLefXGlobal, BottomRightXGlobal)    # Select random X
        random_Y = random.uniform(BottomRightYGlobal, TopLefYGlobal)    # Selext random Y
        randomPnt = [randomID, random_X, random_Y]                      # Points

        distanceFromMiddleOfExtent = distanceCalc2dKm((centerExtentX, centerExtentY), (random_X, random_Y)) 
        
        if distanceFromMiddleOfExtent < r:
            # Check if minimum distance criteria is fulfilled
            if len(randomPnts) > 0:
                dCrit = True
                for i in randomPnts:
                    #distance = math.hypot(random_X - i[0], random_Y - i[1])
                    distance = distanceCalc2dKm([random_X, random_Y], [i[1], i[2]]) # in km
                    if distance <= mindistBetweenWWTPs:  
                        dCrit = False
                if dCrit == True:
                    randomPnts.append(randomPnt)
                    randomID += 1
            else:
                randomPnts.append(randomPnt)
        else:
            _ = 0 #print("too far away")
    return randomPnts, circleExtent

def reDistributeWWTPs(oldWWTP, extent, minDistBetwenWWTP, r):
    '''
    Input: 
    extent                             -    [[topLeftPoint], [topRightPoint]]
    minInterWWTPDistance              -    Minimum Distance Criteria for distribution. The clsoest two point are allowed in initial configuration
    '''
    #print("start redistribution")
    newWWTP = []
    TopLefXGlobal, TopLefYGlobal, BottomRightXGlobal, BottomRightYGlobal = extent[0][0], extent[1][1], extent[1][0], extent[0][1]  # Get Coordinates of Extent
    centerExtentX = TopLefXGlobal + ((BottomRightXGlobal - TopLefXGlobal) / 2)
    centerExtentY = BottomRightYGlobal + ((TopLefYGlobal - BottomRightYGlobal) / 2)

    # Add Distance in order that circle
    TopLefXGlobal = centerExtentX - r * 1000  # in [m]
    TopLefYGlobal = centerExtentY + r * 1000
    BottomRightXGlobal = centerExtentX + r * 1000 
    BottomRightYGlobal = centerExtentY - r * 1000 
    
    circleExtent = [(TopLefXGlobal, BottomRightYGlobal), (BottomRightXGlobal, TopLefYGlobal)]
      
    # Create Random newWWTP
    pos = 0
    while len(newWWTP) < len(oldWWTP):
        random_X = random.uniform(TopLefXGlobal, BottomRightXGlobal)  # Select random X
        random_Y = random.uniform(BottomRightYGlobal, TopLefYGlobal)  # Selext random Y
        distanceFromMiddleOfExtent = distanceCalc2dKm((centerExtentX, centerExtentY), (random_X, random_Y)) 
        
        if distanceFromMiddleOfExtent < r:
            
            # Check if minimum distance criteria is fulfilled
            if len(newWWTP) > 0:
                dCrit = True
                for i in newWWTP:
                    distance = distanceCalc2dKm([random_X, random_Y], [i[1], i[2]]) # in [km]
                    if distance <= minDistBetwenWWTP: 
                        dCrit = False

                if dCrit == True:
                    oldEntry = oldWWTP[pos]
                    oldEntry[1] = random_X
                    oldEntry[2] = random_Y
                    newWWTP.append(oldEntry)
                    pos += 1
            else:
                oldEntry = oldWWTP[pos]
                oldEntry[1] = random_X
                oldEntry[2] = random_Y
                newWWTP.append(oldEntry)
                pos += 1
    return newWWTP, circleExtent

def clusterNN(WWTP, kreisFlache, extent, minInterWWTPDistance, NNATarget, WindowToPlaceRandomNode, NeighbourHoodDensityDistance, NrOfPointsToSelectAndPlaceNewNode):
    '''
    Input: 
    extent                             -    Extent of Circle
    minInterWWTPDistance              -    Minimum Distance Criteria for distribution. The clsoest two point are allowed in initial configuration
    dKritImprov                        -    Neighbourhood for improving index. How far away from node the random node can be selected
    NeighbourHoodDensityDistance       -    How the Density is measured for determinig random point in aggregation
    anzKritNeig                        -    Which density cannot be choosen for node aggregation
    NrOfPointsToSelectAndPlaceNewNode    -    How many nodes are placed in inner window to test
    '''
    
    checked = []
    TopLefXGlobal, TopLefYGlobal, BottomRightXGlobal, BottomRightYGlobal = extent[0][0], extent[1][1], extent[1][0], extent[0][1]  # Get Coordinates of Extent
    randID = len(WWTP) + 1000

    # ------------------------------------------------------------------------------------------------------------------
    # After having created random distributed list with coordinates, iterate until NNA is fulfilled,
    # ------------------------------------------------------------------------------------------------------------------
    #print("INITIAL Circle Area: " + str(kreisFlache))
    clusterDegree = averageNearestNeighborClustering(WWTP, kreisFlache)  # Calculate NNA
    #print("CClusterDegree: " + str(clusterDegree))
    #print("NNATarget: " + str(NNATarget))
    
    if clusterDegree < NNATarget:
        print("The clustering degree is already smaller than the targeted nearest neighbour")
        #print("clusterDegree: " + str(clusterDegree))
        #print("NNATarget: " + str(NNATarget))
        return [] # WWTP
    else:
        listWithNodesTooDense = []  
        while round(clusterDegree, 2) > NNATarget:
            #print("Current NNA-Value: " + str(clusterDegree))
            #print("Nr of checked:     " + str(len(checked)))
            #print("Dense pnts:        " + str(len(listWithNodesTooDense)))
            #print("----------------------------------------------")

            nextPnt, returnToStart = False, 0
            PntsCopy = list(WWTP)
        
            # Find furthest away node
            globalAvarageNNtoAll = 0
            longestSum = 0  # Nearest Distance
            for b in PntsCopy:
                localSum = 0
                for e in PntsCopy:  # Search distance to closest point
                    dst = distanceCalc2dKm((b[1], b[2]), (e[1], e[2]))
                    localSum += dst
                globalAvarageNNtoAll += localSum
                
                if localSum > longestSum and b[0] not in checked:
                    longestSum = localSum
                    XID_fartherstAway = b
            

            if len(checked) == WWTP:
                #print("ERROR: Could not find a distribution")
                return []

            if returnToStart == 0:     
                removedPntList = removePnt(PntsCopy, XID_fartherstAway[0])  # Remove Point from current pnts

                # -------------------------------
                # Get random pnt based on density   
                # -------------------------------
                tooClose = True
                probabilityList = getNrOfNeighboursClustering(PntsCopy, NeighbourHoodDensityDistance)  # Sort list according to number of neighbours 

                # Select random Point and define close range for selection new pnt
                while nextPnt == False and returnToStart == 0: 
                    while tooClose == True and returnToStart == 0:
                       
                        probForDensitSelection = dict(probabilityList)
                        if len(probForDensitSelection) == 0:
                            print("ERROR: densitiyNeighbourhoodToConsiderPacement is too high. No point exists with so many neighbours")
                            checked.append(XID_fartherstAway[0])   
                            returnToStart = 1
                            break

                        # Select node with highest density
                        highestDensit = 0
                        for z in probForDensitSelection:
                            if probForDensitSelection[z] > highestDensit:
                                ranID = z
                                highestDensit = probForDensitSelection[z]
                        
                        # Get pnt coordinate
                        for i in PntsCopy:
                            if i[0] == ranID:
                                windowPnt = i
                                break
 
                        # Remove from probability list
                        scrapCopyList = {}
                        for i in probabilityList:
                            if i != ranID and i not in listWithNodesTooDense:
                                scrapCopyList[i] = probabilityList[i]
                        probabilityList = dict(scrapCopyList)

                        # ---------------------------
                        # Select random pnt in window
                        # ---------------------------
                        #print("WindowToPlaceRandomNode: " + str(WindowToPlaceRandomNode))
                        randomNodeSelection = 0
                        TopLefX = windowPnt[1] - WindowToPlaceRandomNode * 1000         # [m]
                        BottomRightX = windowPnt[1] + WindowToPlaceRandomNode * 1000    # [m]
                        BottomRightY = windowPnt[2] - WindowToPlaceRandomNode * 1000    # [m]
                        TopLefY = windowPnt[2] + WindowToPlaceRandomNode * 1000         # [m]

                        while randomNodeSelection < NrOfPointsToSelectAndPlaceNewNode and returnToStart != 1:  # and selectRandomNode == True:
                            random_X, random_Y = random.uniform(TopLefX, BottomRightX), random.uniform(BottomRightY, TopLefY)    
                            
                            randID = XID_fartherstAway[0]   #
                            randomNEWPnt = XID_fartherstAway
                            randomNEWPnt[0] = randID
                            randomNEWPnt[1] = random_X
                            randomNEWPnt[2] = random_Y
                            randID += 1
                
                            if random_X < TopLefXGlobal or random_X > BottomRightXGlobal or random_Y < BottomRightYGlobal or random_Y > TopLefYGlobal:
                                #print("Out of Extent")
                                _ = 0
                            else: 
                                # -------------------------------------
                                # Check that random node is not too close of any node
                                # -------------------------------------
                                isCloseToBuilding = checkIfIsToocloseToBuilding(removedPntList, minInterWWTPDistance, random_X, random_Y)
    
                                if isCloseToBuilding == True:  # If too close
                                    randomNodeSelection += 1    
                                else:
                                    scrapCopy = list(removedPntList)
                                    scrapCopy.append(randomNEWPnt)
                                    
                                    # takes a lot of time
                                    ANNDistanceNew = averageNearestNeighborClustering(scrapCopy, kreisFlache)
                                    
                                    localGlobalNNAtoAll = 0
                                    longestSum = 0                                 
                                    for b in scrapCopy:
                                        localSum = 0
                                        for e in scrapCopy:                                      # Search distance to closest point
                                            dst = distanceCalc2dKm((b[0], b[1]), (e[0], e[1]))
                                            localSum += dst
                                        localGlobalNNAtoAll += localSum

                                    # check if nearest distance to closer
                                    if ANNDistanceNew < clusterDegree or globalAvarageNNtoAll > localGlobalNNAtoAll:
                                        removedPntList.append(randomNEWPnt)
                                        WWTP = list(removedPntList)
                                        clusterDegree = ANNDistanceNew
                                        nextPnt = True
                                        checked = []
                                        returnToStart = 1 
                                    else:          
                                        print("IMPROVEMENT IS NOT FOUND A")
                                        randomNodeSelection += 1    
                                
                                        # Make faster and test which nodes are already too densifiey
                                        _scrap = NrOfPointsToSelectAndPlaceNewNode  # - 1

                                        # Start kicking in later in algorithm because in advance with innerhouse distance it is soon too dense
                                        if randomNodeSelection == _scrap and returnToStart != 1 and ranID not in listWithNodesTooDense:
                                            listWithNodesTooDense.append(ranID)
                                            returnToStart = 1  
        #print("finished clustering: " + str(clusterDegree))
        return WWTP

def getNrOfNeighbours(PntsCopy, bufferNeighbour):
    probabilityList = {}
    for b in PntsCopy:
        nrOfNeighborsinRadius = 0
        for e in PntsCopy:  # Search distance to closest point
            dst = distanceCalc2dKm((b[0], b[1]), (e[0], e[1]))
            if dst <= bufferNeighbour:
                nrOfNeighborsinRadius += 1
            
        probabilityList[b[0]] = nrOfNeighborsinRadius
    return probabilityList
                
def getNrOfNeighboursClustering(PntsCopy, bufferNeighbour):
    probabilityList = {}
    for b in PntsCopy:
        nrOfNeighborsinRadius = 0
        for e in PntsCopy:  # Search distance to closest point
            dst = distanceCalc2dKm((b[1], b[2]), (e[1], e[2]))
            if dst <= bufferNeighbour:
                nrOfNeighborsinRadius += 1
            
        probabilityList[b[0]] = nrOfNeighborsinRadius
    return probabilityList
                               

def getNearestWWTP(WWTP, depot):
    distInitial = 99999999999
    for i in WWTP:
        currDist = distanceCalc2dKm(depot, (i[0], i[1]))
        if currDist < distInitial:
            distInitial = currDist
            closestWWTP = (i[0], i[1])
    # Del WWTP
    copyWWTP = []
    for i in WWTP:
        if i[0] != closestWWTP[0] and i[1] != closestWWTP[1]:
            copyWWTP.append(i)
    WWTP = copyWWTP
    return WWTP, closestWWTP

def calculateFailureProbability(linareFailureRatePerYear, daysSinceControl, randomAge, dayYear):
    
    failureRateYear = linareFailureRatePerYear
    # failuerRateDay = 1-0.999707 #  0.9=x^360
    failuerRateDay = 1 - math.e ** (math.log(failureRateYear) / dayYear)

    # Straith forward
    failuerRateDay = linareFailureRatePerYear / dayYear
    return failureRateYear, failuerRateDay
  
def addDailyFlow(OST, sludgeAndScumPerPE):
    ''' This function adds the daily flow of a wwtp by multiplying the dimension of a wwtp with average production'''
    for i in OST:
        i[8] += sludgeAndScumPerPE * i[4]       
    return OST

def dimensionWWTP(OST, minWWTPSize, maxWWTPSize, tankSizeINYears, sludgeAndScumPerPE, dayYear, NrYears, givenDistribution, FullFill):
    '''This function assigns a random size and filling status'''

    wwtpID = 0
    newWWTP = []
    for i in OST:
        wwtpID += 1                                             # ID
        Xcor, Ycor, randomZ = i[1], i[2], 100

        if givenDistribution == True:
            randomSize = i[4]
        else:
            randomSize = random.randint(minWWTPSize, maxWWTPSize)  # [EW]
            
        # Tank Size is max for one year 
        randomAge = random.randint(0, 30)  # [year]
        #daysSinceControl = random.randint(0, NrYears * dayYear)  # [day]
        daysLastSmartRoute = random.randint(0, dayYear)  # [day] # Once a year
        
        if FullFill == False:
            randomFill = random.uniform(0, randomSize)  # [EW]

            randomTankFilling = dayYear * tankSizeINYears * sludgeAndScumPerPE * randomFill
            # X, Y, Z, Size, Age, Days since control, working probability, working/YesOrNO, current filling status, days since last smartroute
            newWWTP.append([wwtpID, Xcor, Ycor, randomZ, randomSize, randomAge, 0, True, randomTankFilling, daysLastSmartRoute])  
            
        else:
            fullTank = dayYear * tankSizeINYears * sludgeAndScumPerPE * randomSize
            z = [wwtpID, Xcor, Ycor, randomZ, randomSize, randomAge, 0, True, fullTank, daysLastSmartRoute]
            newWWTP.append(z)
    return newWWTP

def calcAverageOSTSize(WWTP):
    '''Calculate average WWTP Size in PE.'''
    avSizeWWWTP, totalPE = 0, 0
    for i in WWTP:
        avSizeWWWTP += i[4]
        totalPE += i[4]
        
    avSizeWWWTP = float(avSizeWWWTP) / len(WWTP)

    return avSizeWWWTP, totalPE

def getTotalSludgeVolume(totalPE, sludgeAndScumPerPE, dayYear):
    '''Calculate average WWTP Size in PE.'''

    totSludgeScumVolume = totalPE * sludgeAndScumPerPE * dayYear    # Total Sludge and Scum volume
    return totSludgeScumVolume

def checkWWTPtoEmpty(OST, dayYear, tankSizeINYears, sludgeAndScumPerPE, tankFillPeopleCall):
    ''' Check if there is a wwtp to empty. Change criteria that tanks needs emptying
    
    Input
    WWTP                    -    List wiht WWTP
    dayYear                 -    Nr of days a year
    tankSizeINYears         -    Tank size in eary
    sludgeAndScumPerPE      -    accumulation rate per PE
    tankFillPeopleCall      -    Time when people call
    
    Output:
    WWTPstoEmpty            -    All WWTPs which need to be emptied
    OST                    -    Updated WWTP

    '''
    OSTCopy = copy.deepcopy(OST)
    
    WWTPstoEmpty = []
    digitsToRound = 3
    
    for i in OSTCopy:
        totalTankCapacity = i[4] * dayYear * tankSizeINYears * sludgeAndScumPerPE  # Calc total tank capacity
        if i[8] >= round((totalTankCapacity * tankFillPeopleCall), digitsToRound):  # If tank is lager than tank filling at moment in time when people call
            crittanksIsFull = True
        else:
            crittanksIsFull = False
        
        i[7] = crittanksIsFull
        
        if crittanksIsFull == True:
            WWTPstoEmpty.append(i)
            
    
    for i in WWTPstoEmpty:
        for f in OST:
            if i[0] == f[0]:
                f[8] = 0
            
    return WWTPstoEmpty, OST

def assignFailures(WWTP, linareFailureRatePerYear, dayPerYear):
    ''' Define which WWTP have a failure'''
    nrFailes = 0
    for i in WWTP:
        _, failureRateDay = calculateFailureProbability(linareFailureRatePerYear, 1, 1, dayPerYear)  # Test if works or not 
        randomNr = random.uniform(0, 1)
        if randomNr >= failureRateDay: 
            failOrWork = True
        else:
            failOrWork = False
            nrFailes += 1
        i[7] = failOrWork

    failedWWTP = getAllFailures(WWTP)                                           # Get all WWTP which failed

    return WWTP, nrFailes, failedWWTP
    '''
    # In case the sludge needs to be pumped in case of failure. Insert the flow at top position
    FlowFailout = 0 
    for i in WWTP:
        if i[7] == False:
            FlowFailout += i[8]
            WWTPstoEmpty.insert(0, i)   # Priorise and put to first place
    '''

def checkIfthereIsSpaceForAnotherWWTP(currentLoad, cltruck, i):
    ''' 
    check if the next wwtp to check woudl have space in truck
    
    WWTP              -    List with OST
    currentLoad       -    Current Load
    cltruck           -    Truck load capacity
    
    '''
    capacityReached = 0
    
    potentiallyNewLoad = currentLoad                # Potential new load
    stillEmptyCapacity = cltruck - currentLoad      # Potential empty space
    potentiallyNewLoad += i[8]                      # Assumed complete tank is loaded
    
    if potentiallyNewLoad >= cltruck:               # IF more on truck
        capacityReached = 1                         # Max Tank is reached
    return capacityReached, stillEmptyCapacity
   
def getClosestSink(currentPosition, sinks, df):
    shortestDist = 999999999
    
    for i in sinks:
        dist = distanceCalc2dKmFactor((currentPosition[0], currentPosition[1]), (i[0], i[1]), df)
        
        if dist < shortestDist:
            shortestDist = dist
            closestDepot = (i[0], i[1])
    
    return closestDepot

def emptyAtSink(loadedOnTrucks, totWorkTime, avTimeToEmptyAtSink):
    ''' Empty Truck at WWTP. Put load on truck, add time needed'''
    loadedOnTrucks = 0
    totWorkTime += avTimeToEmptyAtSink  # Could be more complex calculation                             
    return loadedOnTrucks, totWorkTime

def getClosestWWTP(WWTP, coordinates):
    ''' This function searches the clostst WWTP and delets it out of the list
    
    Input:    
    WWTP                -    List with WWTP
    coordinates         -    Coordinate to which the closeset WWTP should be found
    
    Output:
    WWTsavingListew     -    List with removed WWTP
    closestWWTP         -    Closest WWTP
    '''
    shortestDist = 999999999999
    WWTsavingListew = []

    for i in WWTP:
        dist = distanceCalc2dKm((i[1], i[2]), (coordinates[1], coordinates[2]))
        if dist < shortestDist:
            shortestDist = dist
            closestWWTP = i

    for i in WWTP:
        if i[1] != closestWWTP[1] and i[2] != closestWWTP[2]:
            WWTsavingListew.append(i)  
  
    return WWTsavingListew, closestWWTP
            
def createRoute(OST, depotCord):
    ''' Search cloest WWTP from Depot, then get next closest until all are visited. Store this order in WWTP
    
    OST    -    list with WWTP
    depotCord   -    Depot
    
    Output:
    visitedWWTP    -    rearranged WWTPs in List according to a nearest neighbour path
    '''
    InitialID = -1
    Depot = (InitialID, depotCord[0], depotCord[1])
    #print("start creating route...")
    if len(OST) == 1:
        return OST
    
    visitedWWTP = []
    
    # Get closest
    while len(OST) >= 1:
        OST, closest = getClosestWWTP(OST, Depot)
        visitedWWTP.append(closest)
        Depot = closest

    return visitedWWTP 

def getAllFailures(WWTP):
    failedWWTP = []
    for i in WWTP:
        if i[7] == False:  # found failure
            failedWWTP.append(i)
    return failedWWTP

def repairWWTP(WWTP, coordinates):
    ''' Repair WWT'''
    for e in WWTP:
        if e[1] == coordinates[0] and e[2] == coordinates[1]:
            e[7] = True  # Repaired
            break          
    return WWTP

def calculateRadius(prozentDistance, tpmax, travelSpeed, maxWWTPSize, cltruck, dayYear, sludgeAndScumPerPE, df):
    ''' calculates Radius and returns distanc in km
    
    Input:
    prozentDistance    -    Percent of max radius
    
    '''
    # timesToEmpty = (float(maxWWTPSize)*dayYear*sludgeAndScumPerPE)/float(cltruck)
    # print("Times To Empty: " + str(timesToEmpty))
    r = ((tpmax / 2) * travelSpeed) / df  # / timesToEmpty  # Because there and back again --> divide by 2
    distanceKm = r * prozentDistance
    return distanceKm

def emptyWWTP(OST, newOST, stillEmptyCapacity, currentLoad, temptyOST, totWorkingTime, currentWorkingHours):
    '''Go to new WWTP and take as much until the vehicle filled.
    
    Input:
    OST                  -    List with WWTP
    newOST                -    New OST coordinates
    stillEmptyCapacity    -    Possible load capacity which still fits the truck
    currentLoad           -    Current load
    temptyOST             -    Time needed to empty truck
    currentWorkingHours   -    How many hours already workign
    '''
    currentWorkingHours += temptyOST
    totWorkingTime += temptyOST
    
    for i in OST:
        if i[1] == newOST[0] and i[2] == newOST[1]:
            currentTank = i[8]
            newFill = currentTank - stillEmptyCapacity
            
            if newFill > 0:                                     # Partial Emptying
                i[8] = newFill
                currentLoad += stillEmptyCapacity   # Truck is full
            else:                                               # Full emptying
                i[8] = 0
                currentLoad += currentTank
    return OST, currentLoad, currentWorkingHours, totWorkingTime
       
def storeOSTDistribution(WWTP):
    ''' This function stores the size-distribution of the WWTPs'''
    OSTDistribution = list(WWTP)
    return OSTDistribution  

def checkIfCurrentWWTPIsEmpty(WWTP, i):
    ''' This function checks if a current OST is empty
    
    Input:
    WWTP              -    List with WWTP
    i                 -    WWTP to check
    
    Output:
    currentIsempty    -    Criteria wheter empty or not
    '''
    for b in WWTP:
        if i[1] == b[1] and i[2] == b[2]:
            if i[8] == 0:
                currentIsEmpty = True
            else:
                #print("CURRENT LOAD: " + str(i[8]))
                currentIsEmpty = False
            return currentIsEmpty
    
def calcYearlyAverage(resultList, yearResults, NrYears, r, circleArea, densityPerKm2, PEperKm2):
    ''' This function calcules the average of many years
    
    resultList     -    List to store results
    yearResults    -    List with results per year
    nrYears        -    Number of years
    r              -    radius
    circleArea     -    circle area
    densityPerKm2  -    Density per km
    PEperKm2       -    Density per PE
    '''

    a, b = [],[]
    #c = []    

    #for iteratio in resultOperationCall:
    for iteratio in yearResults:
        a.append(iteratio[0])
        b.append(iteratio[1])
        #c.append(iteratio[2])   # Scrap
    costs = sum(a) / float(NrYears)
    time = sum(b) / float(NrYears)
            
    z = [r, circleArea, densityPerKm2, PEperKm2, costs, time, 0]    # Average cost of all iterations
    resultList.append(z)
    return resultList

def shrinkRadiusFunction(shrinkRadius, r, currIteration, IterationsToAchieve, Rmin, Rshrink):
    '''
    This function shrinks the radius.
    
    Input:
    shrinkRadius            -    True: Shrink takes place
    r                       -    current radius
    currIteration           -    Current Iteration
    IterationsToAchieve     -    Iterations to achieve
    Rmin                    -    Minimum Radius
    Rshrink                 -    How much the radius gets shrinked
    
    Output:
    r                       -    New Radius
    currIteration           -    Current Iteration
    '''
    if shrinkRadius == True:
        if r == Rmin and currIteration == IterationsToAchieve:  # If current radius is minimum radius and number of iterations achieved exit definitely
            r, currIteration = 0, 99999999999
        else:
            
            if currIteration == IterationsToAchieve:            # If number of iterations achieved, shrink radius
                r -= Rshrink                                    # Decrease Radius
                currIteration = 999999999999999
            if r < Rmin:                                        # Set to minium radius
                r = Rmin                                        # define as minimum
                currIteration = 999999999999999
    return r, currIteration

def convertToCosts(totTravelDistance, totTravelTime, totWorkingTime, laborCostSludgeH, fixedCostKM, crentH, totalPE):
    ''' Convert to Costs 
    
    Input:
    totTravelDistance            -    Tot travel distance
    totTravelTime                -    Tot travel time
    totWorkingTime               -    Tot labor time
    laborCostSludgeH             -    Labor costs per h
    fixedCostKM    -             -    Fixed costs per km
    crentH
    totalPE
    
    Output:
    
    nrOfWeeks, totTimePE, fixedCostOperationPe, variableCostOperationPe, fixedCostOperation, variableCostOperation, sr_tot_operation_costPe         
    
    '''
    # Costs Personell
    costPersonellTravelTime = totTravelTime * laborCostSludgeH              # Cost of driving time for personnel
    costWorkingTime = totWorkingTime * laborCostSludgeH                     # Cost of of working time personnel 
    
    # Cost Vehicle
    costTruckDriveTime = totTravelTime * crentH                             # Variable Costs Truck Renting  
    costTruckRent = totWorkingTime * crentH                                 # Fixed Cost Truck Renting (really without driving time?
    totalCostDistanceKMSludge = fixedCostKM * totTravelDistance             # Travel Cost (distance method)
                
    #print("Anzaghl Stunden gemietet ROUTE: " + str(totWorkingTime + totTravelTime))
    #print("Kosten  Stunden gemietet ROUTE: " + str(costTruckRent + costTruckDriveTime))
    
    # Sum Costs
    totCosts = costPersonellTravelTime + totalCostDistanceKMSludge + costTruckRent + costTruckDriveTime + costWorkingTime
         
    #fixedCostOperation = costTruckRent + costWorkingTime
    #variableCostOperation = float(costPersonellTravelTime + totalCostDistanceKMSludge + costTruckDriveTime)   
    #fixedCostOperationPe = fixedCostOperation / float(totalPE)               # Per PE
    #variableCostOperationPe = variableCostOperation / float(totalPE)         # Per PE
    #print("     ")
    #print("COSTS")
    #print(" Distance Cost: " + str(totalCostDistanceKMSludge))
    #print(" costPersonellTravelTime Cost: " + str(costPersonellTravelTime))
    #print(" costTruckRent Cost: " + str(costTruckRent))
    
    totCostsPE =  totCosts / float(totalPE)        # Per PE
    totTimePE = (totTravelTime + totWorkingTime) / totalPE
    return totTimePE, totCostsPE                        

def convertToCostsserviceTour(totTravelDistance, totTravelTime, totWorkingTime, cp, costPerKM, truckRentPerh, totalPE):
    ''' Convert to Costs 
    
    Input:
    totTravelDistance    -  Total travel distance  
    totTravelTime        -  Total travel time
    totWorkingTime       -  Total working time
    cp                   -  Labor costs 
    costPerKM            -  costs per km 
    truckRentPerh        -  truck rent per h
    totalPE              -  total PE
    
    Output:
    totTimePE        -    Total Time per PE
    totCostPerPE     -    Total Costs per PE
    
    '''
    # Cost Personnel
    costTravelTime = totTravelTime * cp                 # Cost of driving time for personnel
    costWorkingTime = totWorkingTime * cp               # Cost of of working time personnel 
    
    # Cost Vehicle
    costFixedKM =  totTravelDistance * costPerKM         # Travel Cost (distance method)
    costTruckRent = totWorkingTime * truckRentPerh       # Fixed Cost Truck Renting
    costTruckDriveTime = totTravelTime * truckRentPerh   # Variable Costs Truck Renting  
    
    # Sum Costs
    totCosts = costTravelTime + costWorkingTime + costFixedKM  + costTruckRent + costTruckDriveTime # Method fixed cost per km

    # Summen Costs 
    totCostPerPE = totCosts / float(totalPE)  
                    
    # Calculate Parameters
    totTimePE = (totTravelTime + totWorkingTime) / totalPE

    return totTimePE, totCostPerPE     

def writeTotxtSensitivity(outListStep_point, name, inList):
    """
    This functions writes to a .txt file.
    
    Input Arguments: 
    outListStep_point    --    Path to folder where .txt file is stored
    name                 --    Name of .txt file
    inList               --    List with entries to write to .txt file    
    """
    outfile = outListStep_point + name + ".txt"
    myDocument = open(outfile, 'w')
    print(myDocument)
    myDocument.write(str(inList))    
    myDocument.close()
    return

def getRadiusCosts(liste, iterations):
    ''' This function calcules the average costs per radius of the Scheduled Evacuation
    
    Input:
    Liste         -    Results
    iterations    -    Number of iterations
    
    Output:
    avRadiusCosts    -    Average costs per radius
    '''
    avRadiusCosts = []
    
    # Operation Cost Smart Route
    pos, a,b,c = 1, 0,0,0
    std_a, std_b, std_c = [], [], []

    for i in liste:
        a += i[4]
        b += i[5]
        c += i[6]
        std_a.append(i[4])
        std_b.append(i[5])
        std_c.append(i[6])
        
        if pos == iterations:
            pos = 0
            
            # Calculate standard deviations
            SD_a = numpy.std(std_a)
            SD_b = numpy.std(std_b)
            SD_c = numpy.std(std_c)
            
            # Radius, Circle Area, Density per km2, PE per km2, total cost, time, averagecost per m3
            z = [
                 i[0], 
                 i[1], 
                 i[2], 
                 i[3], 
                 a/float(iterations), 
                 b/float(iterations), 
                 c/float(iterations), 
                 SD_a, 
                 SD_b, 
                 SD_c
                ]
            
            avRadiusCosts.append(z)
            a,b,c = 0,0,0
            std_a, std_b, std_c = [], [], []
        pos += 1
    return avRadiusCosts

def summenDayCostsOfYear(costYearCall, dailyCostListCallSludge):
    ''' Summen costs of all days of one year and empty list
    
    Input:
    costYearCall                -    List to store yearly average results
    dailyCostListCallSludge     -    List with costs of each day of one year 
    
    Output:
    costYearCall                -    Update liste with another result of one year
    dailyCostListCallSludge     -    Empty list of one year
    
    '''
    sumTotCosts, sumTotTime = 0, 0
    
    for i in dailyCostListCallSludge:
        sumTotCosts += i[1]             # Cost
        sumTotTime += i[2]              # Time
    
    z = (sumTotCosts, sumTotTime) #, sumTotCostsVariable)
    costYearCall.append(z)
    dailyCostListCallSludge = []                            # Empty list to store daily results           
    return costYearCall, dailyCostListCallSludge

def readInstreetVertices(pathstreetVertices):
    """
    This functions reads in a .txt file.
    
    Input Arguments: 
    pathstreetVertices    --    Path to .txt file
    
    Output Arguments:
    StreetVertices        --    Street Vertices
    """
    txtFileList = readLines(pathstreetVertices)
    StreetVertices = []
    for i in txtFileList:
        lineElements = i.split()
        StreetVertices.append([int(lineElements[0][1:-1]), float(lineElements[1][:-1]), float(lineElements[2][:-1]), float(lineElements[3][0:-1])])
    return StreetVertices

def readInDictionary(pathInFile):
    """
    This functions reads a .txt file into a dictionary.
    
    Input Arguments: 
    pathInFile           --    Path to .txt file
    
    Output Arguments:
    outDictionary        --    Dictionary
    """
    txtFileList = readLines(pathInFile)
    outDictionary = {}  # empty Dictionary
    
    for i in txtFileList:
        spl = i.split(None, 1)  # get index
        index = int(spl[0])
        entries = spl[1].split(",",)
        subDict = {}

        # First Entry
        firstEntry = entries[0].split()
        subDict[int(firstEntry[0][1:-1])] = float(firstEntry[1][:-1])

        # entries in between
        if len(entries) >= 2:
            for entry in entries[1:-1]:
                splitEntry = entry.split(None,)
                subDict[int(splitEntry[0][:-1])] = float(splitEntry[1])
            # Last Entry
            lastEntry = entries[-1].split()
            subDict[int(lastEntry[0][:-1])] = float(lastEntry[1][:-1])
        else:  
            lastEntry = entries[-1].split()
            subDict[int(lastEntry[0][1:-1])] = float(lastEntry[1][:-1])
        outDictionary[index] = subDict
    return outDictionary

def readLines(pathInFile):
    """
    This functions reads out lines of a .txt file
    
    Input Arguments: 
    pathInFile       --    Path to .txt file
    
    Output Arguments:
    readLines        --      Statistics
    """
    inputfile = open(pathInFile, 'r')   # Set Path to existing .txt-File with R results
    lineArray = inputfile.readlines()   # Read in single result lines from txtfile into List
    readLines = []                      # Create empty list to fill in whole txt File with results
    position = 0                        
    while position < len(lineArray):    # Iterate through list with line-results by position
        entry = lineArray[position]     # Get line at position
        readLines.append(entry)         # Append line at position to empty List
        position += 1                   # For the loop
    inputfile.close()                   # Close result .txt file 
    return readLines

def removePnt(PntsCopy, ID):
    '''remove from list'''
    removedPntList = []  
    for e in PntsCopy:
        if e[0] == ID:
            _ = 0
        else:
            removedPntList.append(e)
    return removedPntList

def averageCostPerR(serviceTour, sheduledEvacuation, unsheduledEvacuation, repairTour, resultsService, resultsScheduledEvacuation, resultsUnscheduledEvacuation, resultsRepair, iterations):
    '''
    This function calculates the average for all iterations.
    
    Input:
    serviceTour                    -    Service criteria
    sheduledEvacuation             -    Scheduled Evacuation criteria
    unsheduledEvacuation           -    Unscheduled Evacuation criteria
    repairTour                     -    Repair Tour
    resultsService                 -    Results service
    resultsScheduledEvacuation     -    Results Scheduled Evacuation
    resultsUnscheduledEvacuation   -    Results unscheduled Evacuation
    resultsRepair                  -    Results repair
    iterations                     -    Numer of iterations
    
    Output:
    
    averaged costs
    
    '''
    avRadiusCostserviceTour, avRadiusCostSheduledEvacuation, avRadiusCostUnsheduledEvacuation, avRadiusCostRepair = [], [], [], []
    
    if serviceTour == True:
        #print("Result: Service Tour")
        #avRadiusCostserviceTour = getRadiusCostsserviceTour(resultsService, iterations)
        avRadiusCostserviceTour = getRadiusCosts(resultsService, iterations)
    if sheduledEvacuation == True:
        #print("Result: Sheduled Evacuation")
        #avRadiusCostSheduledEvacuation = getRadiusCostsSRsludge(resultsScheduledEvacuation, iterations)
        avRadiusCostSheduledEvacuation = getRadiusCosts(resultsScheduledEvacuation, iterations)
   
    if unsheduledEvacuation == True:
        #print("RESULT: Unsheduled Evacuation")
        #avRadiusCostUnsheduledEvacuation = getRadiusCostCallSludge(resultsUnscheduledEvacuation, iterations)  
        avRadiusCostUnsheduledEvacuation = getRadiusCosts(resultsUnscheduledEvacuation, iterations)  

    if repairTour == True:
        #print("Result:Repair")
        #avRadiusCostRepair = getRadiusCostCallSludge(resultsRepair, iterations)
        avRadiusCostRepair = getRadiusCosts(resultsRepair, iterations)
        
    # Create a list with summed costs to calculate standard devition of total costs...
    
    return avRadiusCostserviceTour, avRadiusCostSheduledEvacuation, avRadiusCostUnsheduledEvacuation, avRadiusCostRepair

def summenAllCostsScatter(fullCostSumming, IterationsToAchieve, sheduledEvacuation, unsheduledEvacuation, serviceTour, repairTour, listWithIterResultSheduledEvacuation, listWithIterResultsUnsheduled, listWithIterResultServiceTour, listWithIterResultsRepair, totalPE, totSludgeScumVolume):
    
    if fullCostSumming == True:
        fullCostList = []
    
        if sheduledEvacuation == True and unsheduledEvacuation == True:
            print("ERROR: DECIDE WHICH ONE TO WRITE SUMMEN")
            prnt("..")
                    
        # Add sheduled evacuation
        if sheduledEvacuation == True and unsheduledEvacuation == False:
            print("sheduledEvacuation")
            if len(fullCostList) == 0:
                for i in listWithIterResultSheduledEvacuation:
                    fullCostList.append(i)  # Add Smart Route Costs
            else:
                pos = 0
                for i in listWithIterResultSheduledEvacuation:
                    fullCostList[pos][4] += i[4]  # Add Sheduled Cost
                    fullCostList[pos][5] += i[5]  # Add Time
                    fullCostList[pos][6] += i[6]  # scrap
                    pos += 1
        fullCostList_A = copy.deepcopy(fullCostList)
        
        # Add unsheduled evacuation
        if unsheduledEvacuation == True and sheduledEvacuation == False:
            print("unsheduledEvacuation")
            if len(fullCostList_A) == 0:
                for i in listWithIterResultsUnsheduled:
                    fullCostList_A.append(i)  # Add Smart Route Costs
            else:
                pos = 0
                for i in listWithIterResultsUnsheduled:
                    fullCostList_A[pos][4] += i[4]  # Add Unsheduld COst
                    fullCostList_A[pos][5] += i[5]  # Add unsheduled time
                    fullCostList_A[pos][6] += i[6]  # scrap
                    pos += 1
        fullCostList_B = copy.deepcopy(fullCostList_A)         
      
        # Add Service      
        if serviceTour == True:
            if len(fullCostList_B) == 0:
                for i in listWithIterResultServiceTour:
                    fullCostList_B.append(i)        # Add Smart Route Costs
            else:
                pos = 0
                for i in listWithIterResultServiceTour:
                    fullCostList_B[pos][4] += i[4]  # Add service cost
                    fullCostList_B[pos][5] += i[5]  # Add service time
                    fullCostList_B[pos][6] += i[6]  # scrap
                    pos += 1
                    
        # Update if new one
        fullCostList_C = copy.deepcopy(fullCostList_B)      
        
        # Add Repair 
        if repairTour == True:
            if len(fullCostList_C) == 0:
                for i in listWithIterResultsRepair:
                    fullCostList_C.append(i)  # Add Smart Route Costs
            else:
                pos = 0
                for i in listWithIterResultsRepair:
                    fullCostList_C[pos][4] += i[4]  # Add Smart Route Costs
                    fullCostList_C[pos][5] += i[5]  # Add Smart Route Costs
                    fullCostList_C[pos][6] += i[6]  # Add Smart Route Costs
                    pos += 1
                    
        # Update if new one
        fullCostList_Final_scatter = copy.deepcopy(fullCostList_C)   
        
        # -------------------
        # Add Costs per m3
        # -------------------
        for i in fullCostList_Final_scatter:
            averageCostperQubikSludge = (i[4] * totalPE) / totSludgeScumVolume    # Averacge Cost per PE Totale PE  / totalsludgevolume
            i[6] = averageCostperQubikSludge
        
        # -------------------------
        # Calculate Measure for Paper
        # -------------------------
        # Calculate Standard Deviation of the paramters for paper
        measures = []
        a,b,c,std_a, std_b, std_c = 0, 0, 0, [], [], []
        pos = 1
        
        for i in fullCostList_Final_scatter:
            #print("i: " + str(i))
            a += i[4]
            b += i[5]
            c += i[6]
            std_a.append(i[4])              # Cost
            std_b.append(i[5])              # Time
            std_c.append(i[6])              # Cost per m3
        
            if pos == IterationsToAchieve:
                #print("ois;" + str(pos))
                pos = 0
                
                SD_a = numpy.std(std_a)     # Calculate standard deviations of total costs
                SD_b = numpy.std(std_b)     # Calculate standard deviations of used time
                SD_c = numpy.std(std_c)     # Calculate standard deviations of costs per m3
                
                z = [
                     i[0], 
                     i[1], 
                     i[2], 
                     i[3], 
                     a/IterationsToAchieve, # Cost
                     b/IterationsToAchieve, # Time
                     c/IterationsToAchieve, # Cost per m3 
                     SD_a,                  # Standard Devition of total costs
                     SD_b,                  # Standard Devition of total costs
                     SD_c
                    ]
                
                measures.append(z)
                a,b,c = 0,0,0
                std_a, std_b, std_c = [], [], []
            pos += 1
    return fullCostList_Final_scatter, measures

def summenTotalCosts(fullCostSumming, sheduledEvacuation, unsheduledEvacuation, serviceTour, repairTour, avRadiusCostSheduledEvacuation, avRadiusCostUnsheduledEvacuation, avRadiusCostserviceTour, avRadiusCostRepair):
    ''' Summen Costs of different tasks for each Radius (summing averages therfore) '''
    
    if fullCostSumming == True: 
        fullCostList = []           # List with final results
    
        if sheduledEvacuation == True and unsheduledEvacuation == True:
            print("ERROR: DECIDE WHICH ONE TO WRITE SUMMEN")
            prnt("...")
                    
        # Add sludge costs smart route
        if sheduledEvacuation == True and unsheduledEvacuation == False:
            #print("SUMM COSTS: -- Sheduled Costs adding: " + str(avRadiusCostSheduledEvacuation))
            #print(" ")
            if len(fullCostList) == 0:
                for i in avRadiusCostSheduledEvacuation: 
                    fullCostList.append(i)  # Add Smart Route Costs
            else:
                pos = 0
                for i in avRadiusCostSheduledEvacuation: 
                    fullCostList[pos][4] += i[4]  # Add Smart Route Costs
                    fullCostList[pos][5] += i[5]  # Add Smart Route Costs
                    fullCostList[pos][6] += i[6]  # Add Smart Route Costs
                    pos += 1
        fullCostList_A = copy.deepcopy(fullCostList)
   
        # Add sludge costs sludge bedarfsentleerung            
        if unsheduledEvacuation == True and sheduledEvacuation == False:
            #print("SUMM COSTS: -- Unsheduled Costs adding: " + str(avRadiusCostUnsheduledEvacuation))
            #print(" ")
            if len(fullCostList_A) == 0:
                for i in avRadiusCostUnsheduledEvacuation:
                    fullCostList_A.append(i)                # Add Smart Route Costs
            else:
                pos = 0
                for i in avRadiusCostUnsheduledEvacuation:
                    fullCostList_A[pos][4] += i[4]  # Add Smart Route Costs
                    fullCostList_A[pos][5] += i[5]  # Add Smart Route Costs
                    fullCostList_A[pos][6] += i[6]  # Add Smart Route Costs
                    pos += 1
        fullCostList_B = copy.deepcopy(fullCostList_A)         
      
                
        if serviceTour == True:
            #print("SUMM COSTS: -- Service Costs adding: " + str(avRadiusCostserviceTour))
            #print(" ")
            if len(fullCostList_B) == 0:
                for i in avRadiusCostserviceTour: 
                    fullCostList_B.append(i)  # Add Smart Route Costs
            else:
                pos = 0
                for i in avRadiusCostserviceTour: 
                    fullCostList_B[pos][4] += i[4]  # Add Smart Route Costs
                    fullCostList_B[pos][5] += i[5]  # Add Smart Route Costs
                    fullCostList_B[pos][6] += i[6]  # Add Smart Route Costs
                    pos += 1
                    
        # Update if new one
        fullCostList_C = copy.deepcopy(fullCostList_B)      
     
        if repairTour == True:
            #print("SUMM COSTS: -- Repair Costs adding: " + str(avRadiusCostRepair))
            #print(" ")
            if len(fullCostList_C) == 0:
                for i in avRadiusCostRepair: 
                    fullCostList_C.append(i)  # Add Smart Route Costs
            else:
                pos = 0
                for i in avRadiusCostRepair: 
                    fullCostList_C[pos][4] += i[4]  # Add Smart Route Costs
                    fullCostList_C[pos][5] += i[5]  # Add Smart Route Costs
                    fullCostList_C[pos][6] += i[6]  # Add Smart Route Costs
                    pos += 1
                    
        # Update if new one
        fullCostList_Final = copy.deepcopy(fullCostList_C)   
    return fullCostList_Final

# Calibration
def dijkstra(streetNetwork, idp0, idp1):
    """
    This function gets the path from the Djikstra list.

    Input Arguments: 
    streetNetwork             --    Distances to all nodes
    idp0                      --    Start node
    idp1                      --    Endnode
    heightDiff                --    Height Difference
    
    Output Arguments:
    archPathList              --    Updates distances to all nodes.
    distStartEnd              --    Distance between the two nodes.
    """                                              
    distances, listDijkstra = dijkstraAlgorithm(streetNetwork, idp0)                # calculate djikstra distances
    scrapPathDjika, distStartEnd = writePath(listDijkstra, distances, idp0, idp1)                 
    archPathList = archPath(scrapPathDjika, distances)                              # the achPathList contains all intermediary, not pouplated nodes of the street network. List gets afterwards appended to P                
    return archPathList, distStartEnd

def relax(StreetNetwork, u, v, D, P):
    """
    Relaxing function of Djikstra Algorithm

    Input Arguments: 
    StreetNetwork    --    StreetNetwork
    u,v,D            --    Djikstra-related varaibles
    P                --    Djikstra-List with all distances to each node
    """
    inf = float('inf')
    d = D.get(u, inf) + StreetNetwork[u][v] # Possible shortcut estimate
    if d < D.get(v, inf):                   # Is it really a shortcut?
        D[v], P[v] = d, u                   # Update estimate and p

def dijkstraAlgorithm(StreetNetwork, s):
    """
    Djikstra Algorithm to find shortest rout on street network.

    Input Arguments: 
    StreetNetwork    --    StreetNetwork
    s                --    Starting node

    Output Arguments:
    D                --    Distances to every node from s
    P                --    Djikstra-List
    """
    D, P, Q, S = {s:0}, {}, [(0, s)], set()  

    while Q:  # All nodes are checked?
        smallesvalue, cnt = 99999999, 0
        for i in Q:
            if i[0] < smallesvalue:
                smallesvalue, u = i[0], i[1]
                smallesPosition = cnt
            cnt += 1
        del Q[smallesPosition]

        if u in S: continue                     # Already visited? Skip it
        S.add(u)                                # We've visited it now
        for v in StreetNetwork[u]:              # Go through all its neighbors
            relax(StreetNetwork, u, v, D, P)    # Relax the out-edges
            Q.append((D[v], v))                 # visited
    return D, P  

def calculateDistanceFactor(WWTPsmartRoute, streetVertices, streetNetwork, Initialdepot):
    '''
    Compare distance along the street and craw
    '''
    factorCrawDirect = []
    cnt = 0

    InitialdepotiWithID = [0, Initialdepot[0], Initialdepot[1]]
    WWTPsmartRoute.insert(0, InitialdepotiWithID)   # start from depot
    WWTPsmartRoute.append(InitialdepotiWithID)     # return to depot
    
    # Initial find depot
    closestDist = 999999999999
    for i in streetVertices:
        dist = distanceCalc2d([Initialdepot[0], Initialdepot[1]], [i[1], i[2]])
        if dist < closestDist:
                closestTo = i[0]
                closestX = i[1]
                closestY = i[2]
                closestDist = dist
    fromNode = closestTo
    fromX = closestX
    fromY = closestY

    for edge in WWTPsmartRoute:
        closestDist = 999999999999
        for i in streetVertices:
            dist = distanceCalc2d([edge[1], edge[2]], [i[1], i[2]])
            if dist < closestDist:
                closestTo = i[0]
                closestX = i[1]
                closestY = i[2]
                closestDist = dist
        toNode = closestTo
        try:
            _, streetDistanceInM = dijkstra(streetNetwork, fromNode, closestTo)                 # Djikstra   
            streetDistanceInKM = streetDistanceInM / 1000
            crawDistanceinM = distanceCalc2d([fromX, fromY], [closestX, closestY])
            crawDistanceinKM = crawDistanceinM / 1000  
                
            if crawDistanceinM != 0:                      
                f_dist = float(streetDistanceInKM) / float(crawDistanceinKM)
                factorCrawDirect.append(f_dist)

        except KeyError:
            _ = 0

        fromNode = closestTo
        fromX = closestX
        fromY = closestY
        cnt += 1
        
    # Delete first and last
    del WWTPsmartRoute[0]
    del WWTPsmartRoute[len(WWTPsmartRoute)-1]

    # calc average
    summe = 0
    for i in factorCrawDirect:
        summe += i
    if summe > 0:
        averageDistanceFactor = summe / len(factorCrawDirect)
    else:
        averageDistanceFactor = 0
    return averageDistanceFactor


def distanceCalc2d(p0, p1):
    """
    This functions calculates the euclidian distance in 2d space.

    Input Arguments: 
    p0, p1                --    Points with coordinates. Form: (ID, X, Y, Z)
    
    Output Arguments:
    distance              --    Distance between the points not taking in account the height difference.
    """
    distance = math.hypot(p0[0] - p1[0], p0[1] - p1[1])         
    return distance

def writePath(dijkstraList, distancesDijka, start, end):
    """
    This function writes out the direct path from list generated with djikstra algorithm. 
    Start node needs to be the same from where the dijkstra list is calculcated.
    
    Input Arguments: 
    dijkstraList            --    All edges
    distancesDijka          --    List with dijkstra distances 
    start                   --    start node
    end                     --    end node
    
    Output Arguments:
    path                    --    Path
    distStartEnd            --    Total distance from start to end
    """
    path, up, distStartEnd = [end], None, distancesDijka[end]
    
    while up != start and start != end: # As long the startin node is not found
        up = dijkstraList[end]
        path.append(up)
        end = up
    return path, distStartEnd

def archPath(dijkstra, distances):
    """
    This functions writes out the edges from the djikstra algorithm of the new path to be added (archPath).
    
    Input Arguments: 
    dijkstra          --    Path of djikstra
    distances         --    Dijkstra-Distances
    
    Output Arguments:
    archPoints        --    List with path
    """
    archPoints, cnter = [], -1
    dijkstra = dijkstra[::-1]  # Invert pathOrigin to get correct flow to origin
    
    # ArchPath: Calculate Distances on the path and the whole path to archPoints
    for i in dijkstra:
        cnter += 1 
        if cnter > 0:
            nextElement = i
            interDist = distances[nextElement] - distances[oldElement]      
            archPoints.append((oldElement, (nextElement, interDist)))  # Form: [(from, (to, dist)), (....)]
        oldElement = i
    return archPoints