import re
from collections import defaultdict
from PhysicsModel import *


class RA1SusyModel(PhysicsModel):

    def __init__(self):
        self.onlySingleMu = True
        self.binList = []
        self.TFdict = defaultdict(dict)


    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        #Signal strength (POI)
        self.modelBuilder.doVar("r[1.,0.,20.]")
        self.modelBuilder.out.var("r").setConstant(True)
        
        #Get the list of bins with no duplicate
        for aBin in self.DC.bins:            
            if "_had" in aBin: self.binList.append(aBin.replace("_had",""))
            if self.onlySingleMu and ("_mumu" in aBin or "_pho" in aBin): self.onlySingleMu = False # abit of hack here...
        
        print "*** self.binList = ",self.binList
            
        #if self.onlySingleMu: print "### Only single muon CR is used"
        
        #Initialise transfer factors
        for aBin in self.binList:

            print "*** Initialising parameters for bin ",aBin
            exp_ttW_had   = self.DC.exp[aBin+"_had"]["ewk_ttW"]
            exp_Zinv_had  = self.DC.exp[aBin+"_had"]["ewk_Zinv"]
            print "    exp_ttW_had = ",exp_ttW_had
            print "    exp_Zinv_had = ",exp_Zinv_had
            exp_all_mu = 0.
            exp_all_mumu = 0.
            exp_all_pho = 0.

            # f_Zinv
            self.modelBuilder.doVar("f_Zinv_"+aBin+"["+str(exp_Zinv_had/(exp_ttW_had+exp_Zinv_had))+",0.0001,1.]")
            #self.modelBuilder.out.var("f_Zinv_"+aBin).setConstant(True)            

            # ewk_had
            self.modelBuilder.doVar("ewk_had_"+aBin+"["+str(exp_Zinv_had+exp_ttW_had)+","+str(0.8*(exp_Zinv_had+exp_ttW_had))+","+str(1.2*(exp_Zinv_had+exp_ttW_had))+"]")            
            self.modelBuilder.out.var("ewk_had_"+aBin).setConstant(True)            

            # r_bkg            
            self.modelBuilder.doVar("r_bkg_"+aBin+"[1.,0.1,10.]")
            #self.modelBuilder.out.var("r_bkg_"+aBin).setConstant(True)            

            print "    f_Zinv = ",self.modelBuilder.out.var("f_Zinv_"+aBin).getVal()
            print "    ewk_had = ",self.modelBuilder.out.var("ewk_had_"+aBin).getVal()

            # FIXME: should we add signal contamination at some point?
            #TF singlemu
            if aBin+"_mu" in self.DC.exp.keys(): 
                if "ewk_ttW" in self.DC.exp[aBin+"_mu"].keys(): 
                    exp_all_mu += self.DC.exp[aBin+"_mu"]["ewk_ttW"]
                if "ewk_Zinv" in self.DC.exp[aBin+"_mu"].keys(): 
                    exp_all_mu += self.DC.exp[aBin+"_mu"]["ewk_Zinv"]                
                if exp_all_mu > 0:  
                    print "    exp_all_mu = ",exp_all_mu
                    if not self.onlySingleMu: self.TFdict[aBin]["mu"] = exp_ttW_had/exp_all_mu
                    else: self.TFdict[aBin]["mu"] = (exp_ttW_had+exp_Zinv_had)/exp_all_mu
                    self.modelBuilder.doVar("TF_mu_"+aBin+"["+str(self.TFdict[aBin]["mu"])+"]")
                    self.modelBuilder.out.var("TF_mu_"+aBin).setConstant(True)
                    print "    TF_mu = ",self.modelBuilder.out.var("TF_mu_"+aBin).getVal()
                else: 
                    print "Something wrong here in bin ",aBin,": looks like we want to include ",aBin+"_mu","CR but we have no count there!"
                    exit(1)

            else: 
                #print "No single muon control region for this bin",aBin
                self.TFdict[aBin]["mu"] = -99.
                self.modelBuilder.doVar("TF_mu_"+aBin+"["+str(self.TFdict[aBin]["mu"])+"]")
                self.modelBuilder.out.var("TF_mu_"+aBin).setConstant(True)

            #TF doublemu
            if aBin+"_mumu" in self.DC.exp.keys(): 
                if "ewk_ttW" in self.DC.exp[aBin+"_mumu"].keys(): exp_all_mumu += self.DC.exp[aBin+"_mumu"]["ewk_ttW"]
                if "ewk_Zinv" in self.DC.exp[aBin+"_mumu"].keys(): exp_all_mumu += self.DC.exp[aBin+"_mumu"]["ewk_Zinv"]                
                if exp_all_mumu > 0:  
                    self.TFdict[aBin]["mumu"] = exp_Zinv_had/exp_all_mumu
                    self.modelBuilder.doVar("TF_mumu_"+aBin+"["+str(self.TFdict[aBin]["mumu"])+"]")
                    self.modelBuilder.out.var("TF_mumu_"+aBin).setConstant(True)
                else: 
                    print "Something wrong here in bin ",aBin,": looks like we want to include ",aBin+"_mumu","CR but we have no count there!"
                    exit(1)

            else: 
                #print "No double muon control region for this bin",aBin
                self.TFdict[aBin]["mumu"] = -99.
                self.modelBuilder.doVar("TF_mumu_"+aBin+"["+str(self.TFdict[aBin]["mumu"])+"]")
                self.modelBuilder.out.var("TF_mumu_"+aBin).setConstant(True)

            #TF single photon
            if aBin+"_pho" in self.DC.exp.keys(): 
                if "ewk_ttW" in self.DC.exp[aBin+"_pho"].keys(): exp_all_pho += self.DC.exp[aBin+"_pho"]["ewk_ttW"]
                if "ewk_Zinv" in self.DC.exp[aBin+"_pho"].keys(): exp_all_pho += self.DC.exp[aBin+"_pho"]["ewk_Zinv"]                
                if exp_all_pho > 0:  
                    self.TFdict[aBin]["pho"] = exp_Zinv_had/exp_all_pho
                    self.modelBuilder.doVar("TF_pho_"+aBin+"["+str(self.TFdict[aBin]["pho"])+"]")
                    self.modelBuilder.out.var("TF_pho_"+aBin).setConstant(True)
                else: 
                    print "Something wrong here in bin ",aBin,": looks like we want to include ",aBin+"_pho","CR but we have no count there!"
                    exit(1)

            else: 
                #print "No single photon control region for this bin",aBin
                self.TFdict[aBin]["pho"] = -99.
                self.modelBuilder.doVar("TF_pho_"+aBin+"["+str(self.TFdict[aBin]["pho"])+"]")
                self.modelBuilder.out.var("TF_pho_"+aBin).setConstant(True)


            # Define expressions for yields in CR
            self.modelBuilder.factory_("expr::CR_mu_"+aBin+"(\"(1/"+str(exp_all_mu)+")*(1/@0)*(1-@1)*@2*@3\",TF_mu_"+aBin+",f_Zinv_"+aBin+",ewk_had_"+aBin+",r_bkg_"+aBin+")")
            if not self.onlySingleMu:
                self.modelBuilder.factory_("expr::CR_mumu_"+aBin+"(\"("+str(exp_all_mumu)+")*(1/@0)*@1*@2*@3\",TF_mumu_"+aBin+",f_Zinv_"+aBin+",ewk_had_"+aBin+",r_bkg_"+aBin+")")
                self.modelBuilder.factory_("expr::CR_pho_"+aBin+"(\"("+str(exp_all_pho)+")*(1/@0)*@1*@2*@3\",TF_pho_"+aBin+",f_Zinv_"+aBin+",ewk_had_"+aBin+",r_bkg_"+aBin+")")        

        # #This is meant for shape analysis and histograms
        # for aBin in self.binList:
        #     if not self.DC.path_to_file(aBin+"_had","ewk_ttW") == "FAKE":      
        #         ##ttW
        #         tf_ttW        = r.TFile(self.DC.path_to_file(aBin+"_had","ewk_ttW"),"READ")
        #         shape_ttW_had = tf_ttW.Get(self.DC.path_to_shape(aBin+"_had","ewk_ttW"))
        #         shape_ttW_mu  = tf_mu.Get(self.DC.path_to_shape(aBin+"_mu","ewk_ttW"))
        #         #Zinv
        #         tf_Zinv         = r.TFile(self.DC.path_to_file(aBin+"_had","ewk_Zinv"),"READ")
        #         shape_Zinv_had = tf_Zinv.Get(self.DC.path_to_shape(aBin+"_had","ewk_Zinv"))
        #         shape_Zinv_mumu = tf_Zinv.Get(self.DC.path_to_shape(aBin+"_mumu","ewk_Zinv"))
        #         shape_Zinv_pho  = tf_Zinv.Get(self.DC.path_to_shape(aBin+"_pho","ewk_Zinv"))
        #         nMHTbin = range(shape_ttW.GetNbinsX()) # binning should be the same for ttW and Z->inv
        #         for i in nMHTbin:
        #             lowE = str(int(shape_ttW.GetBinLowEdge(i+1)))
        #             #SingleMu TF
        #             exp_ttW_had = shape_ttW_had.GetBinContent(i+1)
        #             exp_all_mu  = shape_ttW_mu.GetBinContent(i+1)+shape_Zinv_mu.GetBinContent(i+1) #should we add also signal?
        #             TFdict[aBin]["mht_"+lowE]["mu"] = exp_ttW_had/exp_all_mu
        #             #DoubleMu TF
        #             exp_Zinv_had  = shape_Zinv_had.GetBinContent(i+1)
        #             exp_all_mumu  = shape_ttW_mumu.GetBinContent(i+1)+shape_Zinv_mumu.GetBinContent(i+1) #should we add also signal?
        #             TFdict[aBin]["mht_"+lowE]["mumu"] = exp_Zinv_had/exp_all_mumu
        #             #SinglePhoton TF
        #             exp_all_pho  = shape_ttW_pho.GetBinContent(i+1)+shape_Zinv_pho.GetBinContent(i+1) #should we add also signal?
        #             TFdict[aBin]["mht_"+lowE]["pho"] = exp_Zinv_had/exp_all_pho

        # print TFdict
        # tf_ttW.Close()
        # tf_Zinv.Close()


        # Define POIs
        poi = 'r'
        self.modelBuilder.doSet("POI",poi)
    
    def getYieldScale(self,bin,process):

        # #DEBUG
        # if self.DC.isSignal[process] == 1: 
        #     #print "Process ",process,"in bin ",bin,"has rate ",self.DC.exp[bin][process]*self.modelBuilder.out.var("r").getVal()
        #     return "r"                     
        # else: return 1
        
        catTag = "_".join(bin.split("_")[:-1])
        regTag = bin.split("_")[-1:][0]
        
        if self.DC.isSignal[process] == 1: 
            print "*** Process ",process,"in bin ",bin,"has rate ",self.DC.exp[bin][process]*self.modelBuilder.out.var("r").getVal()
            return "r"                     
        elif "_had" in bin: 
            print "*** Process ",process,"in bin ",bin,"has rate ",self.DC.exp[bin][process]*self.modelBuilder.out.var("r_bkg_"+catTag).getVal()
            return "r_bkg_"+catTag
        elif "_mu" in bin or "_mumu" in bin or "_pho" in bin: 
            print "*** Process ",process,"in bin ",bin,"has rate ",self.DC.exp[bin][process]*(1/self.modelBuilder.out.var("TF_"+regTag+"_"+catTag).getVal())*(1-self.modelBuilder.out.var("f_Zinv_"+catTag).getVal())*self.modelBuilder.out.var("ewk_had_"+catTag).getVal()
            return "CR_"+regTag+"_"+catTag
        else: 
            print "No rule to scale process ",process,"in bin ",bin
            exit(1)
        
                    


class RA1SusyModel3(PhysicsModel):
    def __init__(self):
        self.jetBins = []
        self.parDict = {}
        self.doShape = False
        self.combinedCards = False
                        
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        # Signal strength (POI)
        self.modelBuilder.doVar("r[1.,0.,20.]")                
        
        # Other parameters (rttW_i, r_Zinv_i)
        for aBin in self.DC.bins:
            if "_had" in aBin: # we do not want to double count the bins
                theIDlist = aBin.split("_")
                if "ht" in theIDlist[0]: #it means we have combined cards
                    self.combinedCards = True
                    if "mht" in theIDlist[3]: self.doShape = True
                    self.jetBins.append(theIDlist[1]+"_"+theIDlist[2])
                else:
                    if "mht" in theIDlist[2]: self.doShape = True
                    self.jetBins.append(theIDlist[0]+"_"+theIDlist[1])
                    
        
        for aCat in self.jetBins:
            if "eq2b" in aCat or "ge3b" in aCat:
                self.modelBuilder.doVar("r_ewk_"+aCat+"[1.,0.1,10.]")
                #self.modelBuilder.out.var("r_ewk_"+aCat).setConstant(True)
                self.parDict[aCat,"ewk_ttW"] = "r_ewk_"+aCat
                self.parDict[aCat,"ewk_Zinv"] = "r_ewk_"+aCat                
            else:
                self.modelBuilder.doVar("r_ewk_ttW_"+aCat+"[1.,0.1,10.]")
                self.modelBuilder.doVar("r_ewk_Zinv_"+aCat+"[1.,0.1,10.]")
                #self.modelBuilder.out.var("r_ewk_ttW_"+aCat).setConstant(True)
                #self.modelBuilder.out.var("r_ewk_Zinv_"+aCat).setConstant(True)
                self.parDict[aCat,"ewk_ttW"] = "r_ewk_ttW_"+aCat
                self.parDict[aCat,"ewk_Zinv"] = "r_ewk_Zinv_"+aCat

        # Define POIs
        poi = 'r'
        self.modelBuilder.doSet("POI",poi)


    def getYieldScale(self,bin,process):

        regTag = bin.split("_")[-1]
        if self.combinedCards: catTag = bin.split("_")[1]+"_"+bin.split("_")[2]
        else: catTag = bin.split("_")[0]+"_"+bin.split("_")[1]
        
        if self.DC.isSignal[process] == 1: return "r"
        else: return self.parDict[catTag,process]


RA1SusyModel = RA1SusyModel()
RA1SusyModel3 = RA1SusyModel3()
