#!/usr/bin/env python
import os, sys, re, optparse,pickle,shutil



if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option("-c","--counter",dest="trigcounter",default='skimAnalyzerCount/SkimReport.pck',help="counter")
    parser.add_option("-s","--sigma",dest="sigma",type=float,help="cross section",default=1.0)

    (options,args) = parser.parse_args()
    #define output dictionary
    output=dict()
    rootFile='vvTreeProducer/tree.root'

    #samples = ['WW','WZ','ZZ']
    samples = ['DYJetsToLL_M50_HT100to200','DYJetsToLL_M50_HT200to400','DYJetsToLL_M50_HT400to600','DYJetsToLL_M50_HT600toInf']

    for sample in samples:
        events = 0
        for directory in [dd for dd in os.listdir('./') if sample+'_Chunk' in dd and not '.pck' in dd and not '.root' in dd]:
            if not os.path.exists(directory+'/'+rootFile):
                continue
            print directory
            if os.path.exists(directory+'/'+options.trigcounter):
                print 'found counter'
                counterFile=open(directory+'/'+options.trigcounter)
                counter=pickle.load(counterFile)
                if len(counter)>1 and counter[1][0]=='Sum Weights':
                    events += counter[1][1]
                else:
                    events += counter[0][1]
            else:
                print 'problem in file/ cannot count events in ',directory
                print 'If is data you can ignore'

        output['events'] = events
        output['weight'] = options.sigma/float(events)
        output['sigma']  = options.sigma

        f=open(sample+".pck","w")
        pickle.dump(output,f)
        f.close()

        #merge root files
        os.system('hadd '+sample+'.root '+sample+'_Chunk*/vvTreeProducer/tree.root')   
     
