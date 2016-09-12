// ------------------------------------------------------------
// This is the main C wrapper for FEWZ
// reads in command line arguments
// calling sequence: ./fewz -i inputfile -h histofile -o outputfile -p datdir
// inputfile: input file containing parameters for the computation
// histofile: input file containing histogramming parameters
// outputfile: results will be printed there (appending LO,NLO,NNLO depending on order)
// datdir: path to dat directory containing grids for PDF sets (it's ok to use . or ..)


#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <time.h>

extern"C"
{
  extern void initf_(char*, int*, char*, int*, char*, int*, int*);
  extern void forall_(int*, int*);
  extern void ending_(char*, int*, char*, int*, double*);
}

/* The following may be useful if one wants to disable checkpointing at 
certain times (I/O?) for Condor */

/* #ifdef CONDORED
extern "C"
{
  extern void ckpt();
  extern void _condor_ckpt_disable();
  extern void _condor_ckpt_enable();
}
#endif */

void usagemessage(char execstr[])
{
	std::cout<<"Usage: "<<execstr<<" -i <input_file> -h <histo_file> -p <pdf_dir> -o <output_file> -l <output_dir> -s <sector number>"<<std::endl;
	std::cout<<"  -i\tinput argument configuration file name, defaulting to \"input.txt\""<<std::endl;
	std::cout<<"  -h\thistogram configuration file name, defaulting to \"histograms.txt\""<<std::endl;
	std::cout<<"  -p\tPDF directory (all data files should be stored in <pdf_dir>/dat/), defaulting to \"..\""<<std::endl;
        std::cout<<"  -o\toutput file name, LO, NLO or NNLO prefix will be added to final file name accordingly, defaulting to \"output.txt\""<<std::endl;
        std::cout<<"  -l\toutput directory (where output file will be written at), defaulting to \".\""<<std::endl;
	std::cout<<"  -s\tSector number (zero indexed)"<<std::endl;
}

int main(int argc, char** argv)
{
  
  char* input;
  char* histo;	
  char* ttmpdir;
  char* output;
  char* outdir;
  int len1,len2,len3,len4,len5,niter,done,isect;
  time_t start,end;
  double dif;

  // Defaults for command line arguments
  output="output.txt";
  input="input.txt";
  histo="histograms.txt";
  ttmpdir="..";
  outdir=".";
  isect=1;

  // Usage prompt if no options
  if ( argc==1 )
    {
      usagemessage(argv[0]);
      exit(0);
    }
  // Process arguments
  for(int k=1; k<argc; ++k)
    {
      if(((argv[k])[0])=='-')
	{
	  if(!strcmp(argv[k],"-i"))
	    {
	      k++;
	      input=argv[k];
	    }
	  else if(!strcmp(argv[k],"-h"))
	    {
	      k++;
              histo=argv[k];
            }
	  else if(!strcmp(argv[k],"-p"))
	    {
	      k++;
	      ttmpdir=argv[k];
	    }
          else if(!strcmp(argv[k],"-o"))
            {
              k++;
              output=argv[k];
            }
          else if(!strcmp(argv[k],"-l"))
            {
              k++;
              outdir=argv[k];
            }
          else if(!strcmp(argv[k],"-s"))
            {
              k++;
              isect=atoi(argv[k])+1; // add 1 because Condor zero-indexed
            }
	  else
	    {
	      std::cout<<"Unknown option: "<<argv[k]<<std::endl;
	      usagemessage(argv[0]);
	      exit(1);
	    }
	}
      else
	{
	  std::cout<<"Unknown option: "<<argv[k]<<std::endl;
          usagemessage(argv[0]);
	  exit(-1);
	}
    }

  (len1)=int(strlen(input));
  (len2)=int(strlen(histo));
  (len3)=int(strlen(ttmpdir));
  (len4)=int(strlen(output));
  (len5)=int(strlen(outdir));

  // Initialization (Fortran)
  initf_(input,&len1,histo,&len2,ttmpdir,&len3,&isect);

      done = 0;
      niter = 0;
      time(&start);
      while(!done)
	{
/*          #ifdef CONDORED
            ckpt();
            _condor_ckpt_disable();
          #endif
*/	
//	Main integration routine (Vegas, will always quit each iteration)  
        forall_(&niter,&done);
/*          #ifdef CONDORED
            _condor_ckpt_enable();
          #endif
*/          niter++;
	  time(&end);
	  dif = difftime (end,start);
//	Write output file (each iteration)
	  ending_(output,&len4,outdir,&len5,&dif);
	}

  return 1;

}

