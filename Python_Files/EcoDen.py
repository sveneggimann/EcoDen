# ======================================================================================
# Copyright 2016  Swiss Federal Institute of Aquatic Science and Technology
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
# Eggimann Sven, Truffer Bernhard, Maurer Max (2016): Economies of density for waste water treatment. Water Research
#
                        
# Contact:   sven.eggimann@eawag.ch
# Version    1.0
# Date:      3.6.2016
# Autor:     Eggimann Sven
# ======================================================================================
import os, sys, random, math, copy, numpy  # Imports

pythonPath = os.getcwd()
sys.path.append(pythonPath)

# Instructions
# -------------
#  1. Change path in line 46, 49, 52-54
#  2. If you want to use other distance weighing factors than in Eggimann et al. 2016, run the script in "calibration mode" by setting calibrateDistanceFactors to True
#     For calibration convert the distance weighing results into cost functions and replace the given functions within the python function getDistanceFactor() 
#  3. Set up the model and decide which costs you want to calculate (lines 145 ff)

calibrateDistanceFactors = False         # False: EcoDen is run with paramters from Eggimann et al. 2016, True: distance calibration mode 

# -----------------------------------------
# Change path according to your installation
# -----------------------------------------

# Path for python files
pythonScriptPath = ".../EcoDenPythonFiles/"    # Path to python files

# Paths used for sensitivity
outPathSensitivity = "../SensitivityResults/"  # Path to store results

# Path to generated .txt files
outListFolder = ".../Folder_containing_street_network/"         # Folder to store Street Network
pathStreetGraph = ".../streetNetwork.txt"                       # Path with stored street .txt from 'street_calibration.py'
pathstreetVertices = ".../streetVertices.txt"                   # Path with stored street .txt from 'street_calibration.py'

sys.path.append(pythonScriptPath)

from functions_EcoDen import *                                                                               # Imports functions

# ----------------------------------------------------------------------------------------------------------------
# Model Parameters
# ----------------------------------------------------------------------------------------------------------------

# Geographic Extension (calculations are assumed on a metric coordinate system).       
maxX, maxY, minX, minY = 500000, 500000, 0, 0                       # [m] The extent given in meters. The Extent must always be larger than circle area. No negative numbers are allowed
depotForStreetNetworkReadIN = [250000, 250000]                      # [Coordinate of  central point of the catchment used for street network calibration (used for transformation)
streetNetworkCenter = [600000, 191000]                              # Actual Coordinates of street network central point (used for transformation)
extent = [[minX, maxY], [maxX, minY]]                               # Extent

# Generic parameters
dayPerYear = 365.0                                                  # [-] Nr of Days of a year

# Model Parameters
y = 5                                                               # [-] Nr of years to iterate for failures
totalNrOfDays = y * dayPerYear                                      # [-] Nr of days in total
tankSizeINYears = 1.0                                               # [year] Tanks Size in years filling
r_acc =  0.25                                                       # [m3/year/PE] Sludge and Scum accumulation rate per year
OsludgePerPEperDay = r_acc / dayPerYear                             # [m3/day/PE] Sludge and Scum accumulation rate per day
p_emp = 0.99                                                        # [%] Probability that a OST is emptied

# Sink Related
temptyTruck = 0.33                                                  # [h] Time to empty truck at sink
pdeposit = 0.015                                                    # [OST/km2]  Sink Density. Note: If smaller than 0.01 takes long time to distribute.

# Personnel Parameters
cp = 7401.0                                                         # [CHF] Monthly salary of driving personnel
cpmonth = cp / (4 * 45.0)                                           # [CHF] Hour wage of personnel
tpmax = 8.0                                                         # [h] Max productive working hours per day    

# Vehicle Parameters
ctruckfix = 33.0/100.00 * 1.87                                      # [CHF/km] # 33l per 100 km * diesel price  
vtruck = 50.0                                                       # [km/h] Travel Speed
ctruckrent = 2250 /  (7.0 * tpmax)                                  # [CHF/h] Cost of renting a sludge truck for a an hour (only working hour). 2250 Euro * PPP factor/ (5 * 8h)
cltruck = 14.0                                                      # [m3] Sludge truck load capacity

# Car Parameters
ccarrent = (4.4*16 + 8*0.8)/24.0                                    # [CHF/h] Cost of truck rent per h (Mobility)
ccarfix = 0.5                                                       # [CHF/km] Cost per km (Mobility) oder 8.0/100.00 * 1.87 --> 8 Liter Durschnittsannahme
vcar = 50                                                           # [km/h] Speed of car

# OST related
pfailure = 0.1                                                      # [%] failure rate per year (probability of functioning plant - failture rate). E.g: 0.99 --> Failure rate of 1%
Odimmax = 20.0                                                      # [PE] Max dimension of OSTTP
Odimmin = 4.0                                                       # [PE] Min dimension of OSTTP
nvisits = 1.0                                                       # [-] Number of visits per year of the OSTs for service
tservice = 0.5                                                      # [h] Average time to service OST
temptyOST = 1.0/3.0                                                 # [h] Needed time to evacuation OST
nOST = 500                                                          # [-] Nr of OST in set P
trepair = 1                                                         # [h] Repair time of OST

# Distribution Parameter for nearest neighbour distribution (they may need to be changed depending on case study)
fNN = 1.0                                                           # [-] Target Distribution of NN Clustering
InitialminInterSinkDistance = 0.05                                  # [km] Distribution parameter
minInterWWTPDistance = 0.001                                        # [km] Distribution parameter
WindowToPlaceRandomNode = 0.25                                      # [km] Distribution parameter
NeighbourHoodDensityDistance = 100.0                                # [km] Distribution parameter
NrOfPointsToSelectAndPlaceNewNode = 4.0                             # [-]  Distribution parameter

# Other
tankFillPeopleCall = 1.0                                            # [%] If nOSTe, only comparable
Rmin = 10.0                                                         # [km] Minimum Radius of Circle
Rshrink = 10.0                                                      # [km] how much r is shrinked iteratively
percentageDistanceFactor = 1                                        # Regular distance factor

# Sensitivity analysis Parameters
sensitivityAnalysis = False                                         # True: Perform Sensitivity Analysis
nrSensitvityRuns = 50                                               # [-] Number of Parameter Configurations for Sensitivity Analysis
sensitivitytpmax = 8.0                                              # [h] Used for initial circle generation for sensitivity analysis
sensitivityvtruck = 50                                              # [km/h] fOR INITIAL Circle generation

# Parameter Ranges or Sensitivity Analysis
pdepositList = [0.015, 0.0075, 0.00325]
temptyTruckList = [0.1666, 0.33, 0.5]
OsludgeList = [0.1/dayPerYear, 0.25/dayPerYear, 0.4/dayPerYear]
vtruckList = [20, 50, 80]
cltruckList = [10, 14, 18]
tpmaxList = [7, 8, 9]
FdList = [0.9, 1.0, 1.1]
nOSTList = [100, 500, 900]
p_empList = [0.99, 0.98, 0.97]
temptyOSTList = [0.1666, 1/3.0, 0.5]
dimensionOSTList = [[4,10], [4,20], [4,30]]        
     
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# Setting up the model - Below the routing algorithm can be selected, whether the model should be run in scheduled or unscheduled mode etc.
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
IterationsToAchieve = 1            # [] Number of times the OST get distributed in the catchment. If you want to calculate standard deviation, must be more than 1

# Evacuation method
scheduledEvacuation = False        # [] True: Calculate Scheduled Evacuation
unscheduledEvacuation = True       # [] True: Calculate Unscheduled Evacuation
    
# Service and repair can be actived
serviceTour = True                 # [] True: Service costs are calculated
repairTour = True                  # [] True: Repair costs are calculated 

# Routing Algorithm
NNAlgorithm = True                 # True: Nearest Neighbour, False: Clarke & Wright

if repairTour == True:
    considerFailedWWTP = True
else: 
    considerFailedWWTP = False
   
# ------------------------------------------------------------------------------------------------------------    
# Read in Street Network for calibration of distance weighting factor
# ------------------------------------------------------------------------------------------------------------
pCnt = 1
sensitvityList, sensitvityListMeasures = [], []                # List to store random parameter configurations for Sensitivity Analysis
    
# Create list with days where average needs to be calculated
daysToWriteOut = []
for i in range(y):
    i += 1
    daysToWriteOut.append(i * dayPerYear) 

# If the distances are to be calibrated
if calibrateDistanceFactors == True: 
    streetNetwork = readInDictionary(pathStreetGraph)                                                                               # Read in street Graph
    streetVerticesWrongCoordinates = readInstreetVertices(pathstreetVertices)                                                       # Read in streetVertices
    streetVertices = transformNetworkToDepot(streetVerticesWrongCoordinates, depotForStreetNetworkReadIN, streetNetworkCenter)      # Transform Network
    print("claibrated read in ")

while pCnt <= nrSensitvityRuns:
    if sensitivityAnalysis == False:
        pCnt = 9999999
    else:
        pCnt += 1
        baseScenario = False
        if baseScenario == False:
        
            # Parameters for Sensitivity analysis
            position_depositDentiy = random.randint(0, 2)
            position_emptyTime = random.randint(0, 2)
            position_WWTPDimension = random.randint(0, 2)
            position_SludgeAndScum = random.randint(0, 2)
            position_TruckSpeed = random.randint(0, 2)
            position_TruckLoad = random.randint(0, 2)
            position_maxWorkH= random.randint(0, 2)
            position_distanceFactor = random.randint(0, 2)
            position_initialOSTSet = random.randint(0, 2)
            position_probabilityatHome = random.randint(0, 2)
            position_temptyOST = random.randint(0, 2)
    
            pdeposit = pdepositList[position_depositDentiy]
            temptyTruck = temptyTruckList[position_emptyTime]
            
            Odimmin = dimensionOSTList[position_WWTPDimension][0]
            Odimmax = dimensionOSTList[position_WWTPDimension][1]
            
            OsludgePerPEperDay = OsludgeList[position_SludgeAndScum]
            vtruck = vtruckList[position_TruckSpeed]                
            cltruck = cltruckList[position_TruckLoad]    
            tpmax = tpmaxList[position_maxWorkH]                           
            nOST = nOSTList[position_initialOSTSet]
            percentageDistanceFactor = FdList[position_distanceFactor]
            p_emp = p_empList[position_probabilityatHome]
            temptyOST = temptyOSTList[position_temptyOST]
        
        print("Anzahl OST to generate. " + str(nOST))
        print("Random Parameters")
        print("Scum and sludge per day:   " + str(OsludgePerPEperDay))
        print("vtruck:                    " + str(vtruck))
        print("cltruck:                   " + str(cltruck))
        print("tpmax:                     " + str(tpmax))
        print("nOST:                      " + str(nOST))
        print("percentageDistanceFactor:  " + str(percentageDistanceFactor))
        print("p_emp:                   " + str(p_emp))
        print("temptyOST:                 " + str(temptyOST))
        print("----------------------------")
        
    # Generate Catchement
    circleRadius = (sensitivityvtruck * (sensitivitytpmax-temptyOST-temptyTruck)) / 2       # Generate the same radius for Sensitivity analysis        
    circleRadius = roundTo10Steps(circleRadius, 10)             # Round radius

    # Limit to maximum size of 180km for sensitivity analysis
    if circleRadius > 180.0:
        circleRadius = 180.0

    OST_DistributionNEU, circleExtent = createRandomDisributionCircle(minInterWWTPDistance, extent, nOST, circleRadius)       # Create Random Distribution in CircleCircle
    OSTDistribution = dimensionWWTP(OST_DistributionNEU, Odimmin, Odimmax, tankSizeINYears, OsludgePerPEperDay, dayPerYear, y, False, True)  # OST are filled to the max
                         
    listWithIterResultsUnsheduled, listWithIterResultsRepair, avRadiusCostUnscheduledEvacuation,  = [], [], []
    listWithDistanceFactors_scheduled, listWithDistanceFactors_unscheduled = [], []

    avSizeWWWTP, totalPE = calcAverageOSTSize(OSTDistribution)                              # Calculate average OST dimension
    totSludgeScumVolume = getTotalSludgeVolume(totalPE, OsludgePerPEperDay, dayPerYear)     # CAlculate Total Sludge and Scum Volume
       
    listWithIterResultscheduledEvacuation = []          # List to store cost of each iteration with full capacity for each iteration [Sludge Collection]
    listWithIterResultServiceTour = []                  # List to store costs of each iteration with full capacity [serviceTour]
    
    # Start calculation
    while circleRadius >= Rmin:
        currIteration = 0                                                                                           # Iteration Initial
        distStreetFactorRouting, distStreetFactorCall = getDistanceFactor(circleRadius, percentageDistanceFactor)   # Get distance factor
        
        print("------------------------------------")
        print("Circle Radius: " + str(circleRadius))
        print("------------------------------------")
        
        while currIteration < IterationsToAchieve:
            day = 0                                     # Day
            dailyCostListUnscheduledSludge = []         # List to store costs of a day of slude collection [Bedarfsentleerung]
            costYear = []                               # List to store average costs of a year
            randomWWTP = []                             # List to store random WWTP
            costYearUnsheduled = []                     # List to store annual costs of call system sludge management
            costYearRepair = []                         # 
            
            call_totTravelDistance, call_tonnenKilometer, call_sr_tot_operation_costPE_call_sum = 0, 0, 0               
            sinks, circleArea = generateSinks(extent, InitialminInterSinkDistance, circleRadius, pdeposit)                 # Generate sinks
            circleExtent = getCircleExtent(extent, circleRadius)

            WWTP = copy.deepcopy(OSTDistribution)                                                                          # Replace WWTP with given Distribution
            avSizeWWWTP, totalPE = calcAverageOSTSize(WWTP)                                                                # Calc average WWTP SIZe
        
            densityPerKm2 = len(WWTP) / circleArea      # OST DEnsity
            PEperKm2 = totalPE / circleArea             # PE Density
            
            print("circleArea:        " + str(circleArea))
            print("densityPerKm2:     " + str(densityPerKm2))
            print("PE Density:        " + str(PEperKm2))
            print("Number of sinks:   " + str(len(sinks)))
            
            # Generate OST Distribution
            randomWWTP = []
            while randomWWTP == []:
                print("Distribute within radius: " + str(circleRadius))
                WWTP = reDistributeWWTPs(WWTP, extent, InitialminInterSinkDistance, circleRadius)                               # Redistribute Dimension of WWTPs
                WWTP = WWTP[0]                                                                        # TODO
                NN_WWTP = clusterNN(WWTP, circleArea, circleExtent, minInterWWTPDistance, fNN, WindowToPlaceRandomNode, NeighbourHoodDensityDistance, NrOfPointsToSelectAndPlaceNewNode)  # HOME D
                
                if NN_WWTP == []:      
                    randomWWTP = []                 # Degree of Clusteirng is not correct
                    print("Repeat distribution as clustering degree is not correct....")
                else:
                    randomWWTP = ["Degree of Clustering is correct"]
                    WWTP = copy.copy(NN_WWTP)
                        
            # Search closest sink and assign Depo
            circleExtentCenter = (((extent[0][1] - extent[0][0]) / 2, (extent[1][0] - extent[1][1])/2))   
            InitialDepot = getClosestSink(circleExtentCenter, sinks, distStreetFactorRouting)

            # --------------------
            # Scheduled Evacuation
            # --------------------
            if scheduledEvacuation == True:
                print("Sludge collection method: Scheduled Evacuation")
                WWTPold = copy.deepcopy(WWTP)
                
                # Create Route
                if NNAlgorithm == True:    
                    scheduledRoute = createRoute(WWTPold, InitialDepot)  
                else:
                    print("start clark global")
                    scheduledRoute = clarkeAndWright(WWTPold, InitialDepot, distStreetFactorRouting)             # Clarke & Wright Route
                    print("end clark global")
                                                
                # Test Distribution and visualize
                #testVisualizeCatchement(scheduledRoute, InitialDepot, sinks)

                # Calculate distance factor (compare street and crow distance)
                if calibrateDistanceFactors == True:
                    print("start calculaing distance factor routing: " + str(len(scheduledRoute)))
                    averageDistanceFactor = calculateDistanceFactor(scheduledRoute, streetVertices, streetNetwork, InitialDepot)
                    listWithDistanceFactors_scheduled.append([circleRadius, "Scheduled Evacuation Tour", averageDistanceFactor])   
                else:
                    listWithDistanceFactors_scheduled.append([circleRadius, "Scheduled Evacuation Tour", 0])    

                # Cost Algorithm Scheduled Evacuation
                totTravelDistance, totTravelTime, totWorkingTime = algorithmEvacuation(scheduledRoute, InitialDepot, tpmax, sinks, vtruck, temptyTruck, distStreetFactorRouting, cltruck, temptyOST, p_emp)
                    
                # Convert Operation Costs into CHF 
                timeSheduled, costSheduledPE = convertToCosts(totTravelDistance, totTravelTime, totWorkingTime, cpmonth, ctruckfix, ctruckrent, totalPE)

                print(" ")
                print("---------------------------------------------")
                print("Costs Sludge Transport [Sheduled Evacuation]")
                print("---------------------------------------------")
                print("totTravelDistance per PE:            " + str(totTravelDistance/totalPE))
                print("tot Travel Distnace:                 " + str(totTravelDistance)) 
                print("totTravelTime:                 [h]:  " + str(totTravelTime))
                print("totWorkingTime:                [h]:  " + str(totWorkingTime))
                print("Total Time                  [h/PE]:  " + str(timeSheduled))
                print("Tot Costs:                [CHF/PE]:  " + str((costSheduledPE)))
                
                # Append to Iteration Full Capacity
                z = [circleRadius,                  # Current number of WWTP
                    circleArea,                     # Circle Area
                    densityPerKm2,                  # Density per km2
                    PEperKm2,                       # PE per km2 (density)
                    costSheduledPE,                 # Tot Scheduled
                    timeSheduled,                   # Tot All Time 
                    0                               # Position for cost per m3
                    ]     
                listWithIterResultscheduledEvacuation.append(z) 
                
            # --------------------
            # Service 
            # --------------------
            if serviceTour == True:        
                print (" ")
                print("Service Tour Calculations...")     
                # In case smart Routing was taken, do not recalculate route
                if scheduledEvacuation == False:
                    aa = copy.deepcopy(WWTP)
                    if NNAlgorithm== True:
                        df = createRoute(aa, InitialDepot)                  
                    else:
                        df = clarkeAndWright(aa, InitialDepot, distStreetFactorRouting)
                    WWTPcontrolRoute = copy.deepcopy(df)
                else:
                    WWTPcontrolRoute = copy.deepcopy(scheduledRoute)

                totTravelDistance, totTravelTime, totWorkingTime = algorithmServiceTour(WWTPcontrolRoute, InitialDepot, InitialDepot, tpmax, vcar, distStreetFactorRouting, tservice)

                # Convert Operation Costs into CHF 
                totTimePerYearPE, totCostserviceTourPE = convertToCostsserviceTour(totTravelDistance, totTravelTime, totWorkingTime, cpmonth, ccarfix, ccarrent, totalPE)

                # Store Costs in List
                z = [circleRadius,                                      # Current nubmer of WWTP
                    circleArea,
                    densityPerKm2,
                    PEperKm2,                                           # PE per km2 (density)
                    totCostserviceTourPE * nvisits,                     # Tot Operation Smart Rout
                    totTimePerYearPE,                                   # Tot fixed costs     # 15.06.2015
                    0                                                   # Empty place for cost per m3
                    ]     
                listWithIterResultServiceTour.append(z)

                '''print("---------------------------------------------")
                print("Costs service Tour")
                print("---------------------------------------------")    
                print("Numer of Visists a year:            " + str(nvisits))
                print("Tot travel Distance          [km]:  " + str(totTravelDistance))        
                print("Tot travel Distance per OST [km/PE]:" + str(totTravelDistance/len(WWTP)))        
                print("Tot Costs:                [CHF/PE]: " + str((totCostserviceTourPE * nvisits)))
                print("Tot Time per PE             [h/PE]: " + str(totTimePerYearPE))
                print("---------------")
                '''

            regulationCostDaysCall, listDistanceFactorsDays = [], [] # List to store results of each day
    
            # used for Unscheduled Evacuation 
            if unscheduledEvacuation == True:
                print("Sludge collection method: Unscheduled Evacuation") 
                WWTP = dimensionWWTP(WWTP, Odimmin, Odimmax, tankSizeINYears, OsludgePerPEperDay, dayPerYear, y, True, False)  # Randomly fill WWTP but take use dimension
   
            failedWWTP, WWTPstoEmpty = [], []  

            # Make copy for unscheduled evacuation
            WWTP_unscheduled = copy.deepcopy(WWTP)
            
            # Iterate over Days
            while day <= totalNrOfDays:
                #print("-----------------DAY Nr: " + str(day) + "   " + str(len(WWTP)))

                WWTP_unscheduled = addDailyFlow(WWTP_unscheduled, OsludgePerPEperDay)  # Update flow. Add daily generated flow to WWTP
                
                # ----------------------------------------------------
                # Repair Costs
                # ----------------------------------------------------
                if considerFailedWWTP == True:
                    WWTP_unscheduled, NrOfFailes, failedWWTP = assignFailures(WWTP_unscheduled, pfailure, dayPerYear)  # Assign failures

                    # Create NN Route
                    if NNAlgorithm == True:
                        failedWWTP = createRoute(failedWWTP, InitialDepot)                      # Smart NN-Routing
                    else:
                        if len(failedWWTP) > 2:
                            failedWWTP = clarkeAndWright(failedWWTP, InitialDepot, distStreetFactorCall)
                    
                # Each day the total distance and time is calculated and then summed over the whole year.
                if repairTour == True and len(failedWWTP) > 0:   
                    # In case smart Routing was taken, do not recalculate route
                    if scheduledEvacuation == False:
                        if NNAlgorithm == True:
                            df = createRoute(failedWWTP, InitialDepot)    
                            WWTPfailureRoute = copy.deepcopy(df)   
                        else:
                            if len(failedWWTP) > 2:
                                WWTPfailureRoute = clarkeAndWright(failedWWTP, InitialDepot, distStreetFactorCall)
                            else:
                                WWTPfailureRoute = failedWWTP
                    else:
                        WWTPfailureRoute = copy.deepcopy(failedWWTP)
    
                    totTravelDistance, totTravelTime, totWorkingTime = algorithmServiceTour(WWTPfailureRoute, InitialDepot, InitialDepot, tpmax, vcar, distStreetFactorCall, trepair)

                    # Convert Operation Costs into CHF                                                                                                                                                              
                    totTimePerYearPE, totCostserviceTourPE = convertToCostsserviceTour(totTravelDistance, totTravelTime, totWorkingTime, cpmonth, ccarfix, ccarrent, totalPE)
        
                    '''print("---------------------------------------------")
                    print("Costs Repairing WWTPs " + str(totalPE))
                    print("---------------------------------------------")    
                    print("Tot travel Distance       [km]: " + str(totTravelDistance))        
                    print("totTravelTime:                : " + str(totTravelTime))
                    print("totWorkingTime:               : " + str(totWorkingTime))
                    print("totCostserviceTourPE:           " + str(totCostserviceTourPE))
                    print("totCostserviceTour:           " + str(totCostserviceTourPE * totalPE))
                    print("---------------")
                    '''
                    
                    # Store Costs in List
                    z = [day,                       # Current nubmer of days
                        totCostserviceTourPE,       # Cost Service
                        totTimePerYearPE,           # Time Service
                        0                           # Position for cost per m3
                        ]
                    regulationCostDaysCall.append(z)               
    
                # --------------------------------------------
                # Unscheduled Evacuation 
                # --------------------------------------------           
                if unscheduledEvacuation == True:
                    #print("Unscheduled Evacuation " + str(len(WWTP_unscheduled)))
                    
                    # Get WWTP to empty and empty in WWTP
                    WWTPstoEmpty, WWTP_unscheduled = checkWWTPtoEmpty(WWTP_unscheduled, dayPerYear, tankSizeINYears, OsludgePerPEperDay, tankFillPeopleCall)  # Test if there is WWTP to empty
                    
                    # Check how many wwtps there are to empty
                    if len(WWTPstoEmpty) == 0: 
                        continueEmpting = 0
                    else:
                        # Create NNAlgorithm Route
                        if NNAlgorithm == True:
                            WWTPstoEmpty = createRoute(WWTPstoEmpty, InitialDepot)  # If more than one WWTP per day, create route
                        else:
                            if len(WWTPstoEmpty) > 2: 
                                WWTPstoEmpty = clarkeAndWright(WWTPstoEmpty, InitialDepot, distStreetFactorCall)
                        
                        # Calculate distance factor (compare street and craw distance)
                        if calibrateDistanceFactors == True:
                            averageDistanceFactorSludge = calculateDistanceFactor(WWTPstoEmpty, streetVertices, streetNetwork, InitialDepot)
                            listDistanceFactorsDays.append([circleRadius, "Scheduled Method", averageDistanceFactorSludge]) 
                        else:
                            listDistanceFactorsDays.append([circleRadius, "Scheduled Method", 0])   

                        # Cost Algorithm unscheduled
                        totTravelDistance, totTravelTime, totWorkingTime = algorithmEvacuation(WWTPstoEmpty, InitialDepot, tpmax, sinks, vtruck, temptyTruck, distStreetFactorCall, cltruck, temptyOST, p_emp)
                        #print("start control f...")

                        controlA(WWTPstoEmpty)              # Control if all WWTPs have been emptied
                        
                        # print("totTravelDistance: " + str(totTravelDistance))
                        call_totTravelDistance += totTravelDistance
    
                        # Convert to CHF
                        totTimePE, totCostScheduledEvacuation = convertToCosts(totTravelDistance, totTravelTime, totWorkingTime, cpmonth, ctruckfix, ctruckrent, totalPE) 
    
                        # Append to Iteration Full Capacity
                        z = [day,                               # Day
                            totCostScheduledEvacuation,         # Cost Sheduled
                            totTimePE                           # Time Sheduled
                            ]
                        dailyCostListUnscheduledSludge.append(z)      
                        call_sr_tot_operation_costPE_call_sum += totCostScheduledEvacuation # SCrap  
                        
                        '''print(" ")
                        print("---------------------------------------------")
                        print("Costs Sludge Transport [DAILY ITERATION]")
                        print("---------------------------------------------")   
                        print("tot Travel Distnace in this route:                 " + str(totTravelDistance))       
                        print("TONNENKILOMETER:                     " + str(tonnenKilometer)) 
                        print("totTravelTime:                    [h]" + str(totTravelTime))
                        print("totWorkingTime:                   [h]" + str(totWorkingTime))
                        print("Tot Costs:                [CHF/PE]:  " + str((totCostScheduledEvacuation)))
                        print("Variable Costs:           [CHF/PE]:  " + str(variableCostOperationPe))
                        print("Fixed Costs               [CHF/PE]:  " + str(fixedCostOperationPe))
                        print("---------------")
                        '''
               
                # Add Yearly Costs     
                if day in daysToWriteOut:  
                    # Summen Unsheduled Costs
                    if unscheduledEvacuation == True: 
                        costYearUnsheduled, dailyCostListUnscheduledSludge = summenDayCostsOfYear(costYearUnsheduled, dailyCostListUnscheduledSludge)
                    # Summen Repair Tour Costs
                    if repairTour == True:
                        costYearRepair, regulationCostDaysCall = summenDayCostsOfYear(costYearRepair, regulationCostDaysCall)
                day += 1  # Add Day
            
            #print("----------------------------------- ----------------------------------------------")
            #print("Total Distance scheduled:                                     " + str(call_totTravelDistance/y))
            #print("Total costs scheduled per year                                " + str(call_sr_tot_operation_costPE_call_sum/y))
            #print("Total KM per year :                                           " + str(call_totTravelDistance/y))
            
            if unscheduledEvacuation == True:
                av = 0
                for i in listDistanceFactorsDays:
                    av += i[2]
                if av > 0:
                    distanceFactorSludgeCalling = av / float(len(listDistanceFactorsDays))
                    listWithDistanceFactors_unscheduled.append([circleRadius, "Unscheduled System", distanceFactorSludgeCalling]) 
            
            # Calculate average of all years
            if unscheduledEvacuation == True:  
                listWithIterResultsUnsheduled = calcYearlyAverage(listWithIterResultsUnsheduled, costYearUnsheduled, y, circleRadius, circleArea, densityPerKm2, PEperKm2)
                
            if repairTour == True:
                listWithIterResultsRepair = calcYearlyAverage(listWithIterResultsRepair, costYearRepair, y, circleRadius, circleArea, densityPerKm2, PEperKm2)
            currIteration += 1
                
            # Check if all iterations are reached and radius needs to be shrinked 
            circleRadius, currIteration = shrinkRadiusFunction(True, circleRadius, currIteration, IterationsToAchieve, Rmin, Rshrink)
               
    # Calculate average for each Radius
    avRadiusCostserviceTour, avRadiusCostscheduledEvacuation, avRadiusCostUnscheduledEvacuation, avRadiusCostRepair = averageCostPerR(serviceTour, scheduledEvacuation, unscheduledEvacuation, repairTour, listWithIterResultServiceTour, listWithIterResultscheduledEvacuation, listWithIterResultsUnsheduled, listWithIterResultsRepair, IterationsToAchieve)
    
    # --------------------
    # Printing out results
    # --------------------
    print(" ")
    print(" =================================================")
    print(" Results ")
    print(" =================================================")
    print(" ")
    if calibrateDistanceFactors == True:
        print(" Distance factors unscheduled routing")
        print(" ------------------------------------ ")
        print(" [Radius, scheduling mode, distance factor]")
        for i in listWithDistanceFactors_unscheduled:
            print(i)
        print(" ")
        print(" Distance factors scheduled routing")
        for i in listWithDistanceFactors_scheduled:
            print(i)
        print(" ")

    print(" ")
    print(" Average costs for each radius for scheduled evacuation ")
    print(" ------------------------------------ ")
    print(" Results in List: [Catchment Radius, Catchment Area, Density per km2, PE per km2, Total costs per PE, Total time per PE, costs per m3, Standard Deviation Total costs per PE, Standard Deviation Total time per PE, _]")
    #print(avRadiusCostscheduledEvacuation)
    for i in avRadiusCostscheduledEvacuation:
        print(i)
    print(" ")
    print(" ")
    print(" Average costs for each radius for unscheduled evacuation ")
    print(" ------------------------------------ ")
    print(" Results in List: [Catchment Radius, Catchment Area, Density per km2, PE per km2, Total costs per PE, Total time per PE, costs per m3, Standard Deviation Total costs per PE, Standard Deviation Total time per PE, _]")
    for i in avRadiusCostUnscheduledEvacuation:
        print(i)
    print(" ")
    print(" ")
    print(" Average costs for each radius for service ")
    print(" ------------------------------------ ")
    print(" Results in List: [Catchment Radius, Catchment Area, Density per km2, PE per km2, Total costs per PE, Total time per PE, costs per m3, Standard Deviation Total costs per PE, Standard Deviation Total time per PE, _]")
    for i in avRadiusCostserviceTour:
        print(i)
    print(" ")
    print(" ")
    print(" Average costs for each radius for repair ")
    print(" ------------------------------------ ")
    print(" Results in List: [Catchment Radius, Catchment Area, Density per km2, PE per km2, Total costs per PE, Total time per PE, costs per m3, Standard Deviation Total costs per PE, Standard Deviation Total time per PE, _]")
    for i in avRadiusCostRepair:
        print(i)
    
    # Summen Costs of different tasks for each Radius (summing averages therfore)
    fullCostList_Final = summenTotalCosts(True, scheduledEvacuation, unscheduledEvacuation, serviceTour, repairTour, avRadiusCostscheduledEvacuation, avRadiusCostUnscheduledEvacuation, avRadiusCostserviceTour, avRadiusCostRepair)
    
    print(" ")
    print(" ")
    print(" Average total costs for each radius")
    print(" ------------------------------------ ")
    print(" Results in List: [Catchment Radius, Catchment Area, Density per km2, PE per km2, Total costs per PE, Total time per PE, costs per m3, Standard Deviation Total costs per PE, Standard Deviation Total time per PE, _]")
    for i in fullCostList_Final:
        print(i)

    # Write out all points, not aggregated to radius but summarized Costs and measures
    fullCostList_Final_scatter, totCostsWithCorrectSD = summenAllCostsScatter(True, IterationsToAchieve, scheduledEvacuation, unscheduledEvacuation, serviceTour, repairTour, listWithIterResultscheduledEvacuation, listWithIterResultsUnsheduled, listWithIterResultServiceTour, listWithIterResultsRepair, totalPE, totSludgeScumVolume)
    
    print(" ")
    print(" ")
    print(" Average Total Costs per radius with standard deviations")
    print(" ------------------------------------ ")
    print(" Results in List: [Catchment Radius, Catchment Area, Density per km2, PE per km2, Total costs per PE, Total time per PE, costs per m3, Standard Deviation Total costs per PE, Standard Deviation Total time per PE, _]")
    for i in totCostsWithCorrectSD:
        print(i)
    
    # Sensitivity Analysis  
    if sensitivityAnalysis == True:
        informationSensitivity = [pdeposit, temptyTruck, nOST, OsludgePerPEperDay * dayPerYear, vtruck, cltruck, tpmax, percentageDistanceFactor, temptyOST, Odimmin, Odimmax] # All space sentivie parameters
        sensitvityList.append([informationSensitivity, fullCostList_Final])             # Store whole result in List
        sensitvityListMeasures.append([informationSensitivity, totCostsWithCorrectSD])  # Correct SD for total costs
        print("len(sensitvityList): " + str(len(sensitvityList)))
        
        writeTotxtSensitivity(outPathSensitivity, "_finalSensitivity_50runs_NN_unscheduled", sensitvityList)                     # Write to file
        writeTotxtSensitivity(outPathSensitivity, "_finalSensitivity_50runs_NN_unscheduled_measures", sensitvityListMeasures)    # Write to file

        print(" Results Sensitivity Analysis")
        print(" ----------------------------")
        print(sensitvityList)
