from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
import PhysicsTools.HeppyCore.framework.config as cfg


class XZZDumpEvtList( Analyzer ):
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(XZZDumpEvtList,self).__init__(cfg_ana, cfg_comp, looperName)


    def declareHandles(self):
        super(XZZDumpEvtList, self).declareHandles()


    def beginLoop(self, setup):
        super(XZZDumpEvtList,self).beginLoop(setup)

        self.nevts = 0

        self.outfile = open(self.dirName+'/eventlist.txt', 'w') 

        #line = '*'+'*'*11*4
        #self.outfile.write(line+'\n')
        #line = '*'+'Row'.rjust(9)+' *'+'run'.rjust(9)+' *'+'lumisec'.rjust(9)+' *'+'event'.rjust(9)+' *'
        #self.outfile.write(line+'\n')
        #line = '*'+'*'*11*4
        #self.outfile.write(line+'\n')


    def process(self, event):
        self.readCollections( event.input )

        #line  = '*'+str(self.nevts).rjust(9)+' *'
        #line += str(event.input.eventAuxiliary().id().run()).rjust(9)+' *'
        #line += str(event.input.eventAuxiliary().id().luminosityBlock()).rjust(9)+' *'
        #line += str(event.input.eventAuxiliary().id().event()).rjust(9)+' *'

        line  = str(event.input.eventAuxiliary().id().run())+','
        line += str(event.input.eventAuxiliary().id().luminosityBlock())+','
        line += str(event.input.eventAuxiliary().id().event())

        self.outfile.write(line+'\n')

        self.nevts += 1        

    def endLoop(self, setup):
        super(XZZDumpEvtList,self).endLoop(setup)

        #line = '*'+'*'*11*4
        #self.outfile.write(line+'\n')
        self.outfile.close()


