import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
import json


class XZZDataCardMaker:
    def __init__(self,finalstate,category,luminosity=1.0,physics="LJ"):
        self.physics=physics
        self.finalstate=finalstate
        self.category=category
        self.contributions=[]
        self.systematics=[]
        self.observation=0.0

        self.tag=self.physics+"_"+finalstate+"_"+category
        self.luminosity=luminosity


    def addSystematic(self,name,kind,values,addPar = ""):
        self.systematics.append({'name':name,'kind':kind,'values':values })

    def addObservation(self,observation=0):
        self.observation=observation

    def addContribution(self,name,ID,rate):
        self.contributions.append({'name':name,'ID':ID,'rate':rate}) 

    def makeCard(self):

        f = open(self.tag+'.txt','w')
        f.write('imax 1\n')
        f.write('jmax {n}\n'.format(n=len(self.contributions)-1))
        f.write('kmax *\n')
        f.write('-------------------------\n')
        f.write('bin '+self.tag+'\n')
        f.write('observation '+str(self.observation)+'\n')
        f.write('-------------------------\n')

        # sorted contributions
        contributions = sorted(self.contributions,key=lambda x: x['ID'])
 
        # print bin
        f.write('bin\t') 
        for contrib in contributions:
            f.write(self.tag+'\t')
        f.write('\n')

        #print names
        f.write('process\t')
        for contrib in contributions:
            f.write(contrib['name']+'\t')
        f.write('\n')
       
        #print IDs
        f.write('process\t')
        for contrib in contributions:
            f.write(str(contrib['ID'])+'\t')
        f.write('\n')
 
        #print rate
        f.write('rate\t')
        for contrib in contributions:
            f.write(str(contrib['rate'])+'\t')
        f.write('\n')

        f.write('-------------------------\n')

        # print systematics
        for syst in self.systematics:
            if syst['kind'] == 'param':
                f.write(syst['name']+'\t'+'param\t' +str(syst['values'][0])+'\t'+str(syst['values'][1])+'\n')
            elif syst['kind'] == 'lnN': 
                f.write(syst['name']+'\t'+ 'lnN\t' )
                for contrib in contributions:
                    has=False
                    for name,v in syst['values'].iteritems():
                        if contrib['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )
                            
                        
        f.close()


        

