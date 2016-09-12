import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
import json


class XZZDataCardMaker:
    def __init__(self,channels,luminosity=1.0,tag="LJ"):
        self.tag=tag
        self.channels=channels
        self.contributions=[]
        self.systematics=[]
        self.observations={}

        self.luminosity=luminosity


    def addSystematic(self,name,kind,values,addPar = ""):
        self.systematics.append({'name':name,'kind':kind,'values':values })

    def addObservation(self,observations):
        self.observations=observations

    def addContribution(self,name,ID,rate):
        self.contributions.append({'name':name,'ID':ID,'rate':rate}) 

    def makeCard(self):

        f = open(self.tag+'.txt','w')
        f.write('imax {n}\n'.format(n=len(self.channels)))
        f.write('jmax {n}\n'.format(n=len(self.contributions)-1))
        f.write('kmax *\n')
        f.write('-------------------------\n')
        f.write('bin\t')
        for channel in self.channels:
            f.write(channel+'\t')
        f.write('\n')
        f.write('observation\t')
        for channel in self.channels:
            f.write(str(self.observations[channel])+'\t')
        f.write('\n')
        f.write('-------------------------\n')

        # sorted contributions
        contributions = sorted(self.contributions,key=lambda x: x['ID'])
 
        # print bin
        f.write('bin\t')
        for ch in self.channels:
            for contrib in contributions:
               if ch in contrib['rate'].keys():
                   f.write(ch+'\t')
        f.write('\n')
        
        # print names
        f.write('process\t')
        for ch in self.channels:
            for contrib in contributions:
                if ch in contrib['rate'].keys():
                    f.write(contrib['name']+'\t')
        f.write('\n')
       
        #print IDs
        f.write('process\t')
        for ch in self.channels:
            for contrib in contributions:
                if ch in contrib['rate'].keys():
                    f.write(str(contrib['ID'])+'\t')
        f.write('\n')
 
        #print rate
        f.write('rate\t')
        for ch in self.channels:
            for contrib in contributions:
                if ch in contrib['rate'].keys():
                    f.write(str(contrib['rate'][ch])+'\t')
        f.write('\n')

        f.write('-------------------------\n')

        # print systematics
        for syst in self.systematics:
            if syst['kind'] == 'param':
                f.write(syst['name']+'\t'+'param\t' +str(syst['values'][0])+'\t'+str(syst['values'][1])+'\n')
            elif syst['kind'] == 'lnN': 
                f.write(syst['name']+'\t'+ 'lnN\t' )
                for channel in self.channels:
                    if channel in syst['values'].keys():
                        for contrib in contributions:
                            if contrib['name'] in syst['values'][channel].keys():
                                f.write(str(syst['values'][channel][contrib['name']])+'\t' )
                            else:
                                f.write('-\t' )
                    else:
                        for contrib in contributions:
                            f.write('-\t' )
                    # end if channel.. else ..
                # end for channel ..
            # end if syst['kind'] .. elif ...
            f.write('\n' )
        # end for syst in ...

        # finish card
        f.close()

## end 


###
#                for channel in self.channels:
#                    hasChal=False
#                    for ch,vals in syst['values'].iteritems():
#                        if channel==ch: hasChal=True
#                    if hasChal:
#                        for contrib in contributions:
#                            has=False
#                            for name,v in syst['values'][channel].iteritems():
#                                if contrib['name']==name:
#                                    f.write(str(v)+'\t' )
#                                    has=True
#                                    break;
#                            if not has:
#                                f.write('-\t' )
#                        # end for contributions
#                    else:
#                        for contrib in contributions:
#                            f.write('-\t' )
#                    #end if hasChal
#                # end for channels 
       
#                f.write('\n' )
#            # end if syst['kind']...                
#                        
#        f.close()


        

