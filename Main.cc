/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "Aiger.h"
#include "graph/GraphParser.h"
#include "core/Dimacs.h"
#include "core/Solver.h"
#include "Aiger.h"
#include "core/Config.h"
#include <unistd.h>
using namespace Minisat;

//=================================================================================================


void printStats(Solver& solver)
{
    double cpu_time = cpuTime();
    double mem_used = memUsedPeak();
    printf("restarts              : %" PRIu64 "\n", solver.starts);
    printf("conflicts             : %-12" PRIu64 "   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
    printf("decisions             : %-12" PRIu64 "   (%4.2f %% random) (%.0f /sec)\n", solver.decisions, (float)solver.rnd_decisions*100 / (float)solver.decisions, solver.decisions   /cpu_time);
    printf("propagations          : %-12" PRIu64 "   (%.0f /sec)\n", solver.propagations, solver.propagations/cpu_time);
    printf("conflict literals     : %-12" PRIu64 "   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    if (mem_used != 0) printf("Memory used           : %.2f MB\n", mem_used);
    printf("CPU time              : %g s\n", cpu_time);
}


static Solver* solver;
// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt();   fflush(stdout);}

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
    printf("\n"); printf("*** INTERRUPTED ***\n");
    if (solver->verbosity > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }
    fflush(stdout);
    _exit(1); }


//=================================================================================================
// Main:


int main(int argc, char** argv)
{
    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        // printf("This is MiniSat 2.0 beta\n");
        
#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        printf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));
        
        IntOption    opt_k("bmc", "k","Maximum number of steps to unroll.\n", INT32_MAX, IntRange(0, INT32_MAX));

        StringOption    opt_graph("graph", "graph","", "");


        parseOptions(argc, argv, true);


        mincutalg = ALG_EDMONSKARP;

        if(!strcasecmp(opt_min_cut_alg,"ibfs")){
        	mincutalg=ALG_IBFS;

        }else if (!strcasecmp(opt_min_cut_alg,"edmondskarp")){
        	mincutalg = ALG_EDMONSKARP;
        }else{
        	fprintf(stderr,"Error: unknown max-flow/min-cut algorithm %s, aborting\n",((string)  opt_min_cut_alg).c_str());
        	exit(1);
        }

        reachalg = ALG_CONNECTIVITY;

		 if(!strcasecmp(opt_reach_alg,"dijkstra")){
			reachalg=ALG_DIJKSTRA;

		 }else if(!strcasecmp(opt_reach_alg,"bfs")){
			reachalg=ALG_BFS;

		 }else if (!strcasecmp(opt_reach_alg,"connectivity")){
			 reachalg = ALG_CONNECTIVITY;
		 }else{
			fprintf(stderr,"Error: unknown reachability algorithm %s, aborting\n", ((string) opt_reach_alg).c_str());
			exit(1);
		 }


        double initial_time = cpuTime();

        const char * graphstr = opt_graph;


        // Use signal handlers that forcibly quit until the solver will be able to respond to
        // interrupts:
        signal(SIGINT, SIGINT_exit);
        signal(SIGXCPU,SIGINT_exit);

        // Set limit on CPU-time:
        if (cpu_lim != INT32_MAX){
            rlimit rl;
            getrlimit(RLIMIT_CPU, &rl);
            if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
                rl.rlim_cur = cpu_lim;
                if (setrlimit(RLIMIT_CPU, &rl) == -1)
                    printf("WARNING! Could not set resource limit: CPU-time.\n");
            } }

        // Set limit on virtual memory:
        if (mem_lim != INT32_MAX){
            rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
            rlimit rl;
            getrlimit(RLIMIT_AS, &rl);
            if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
                rl.rlim_cur = new_mem_lim;
                if (setrlimit(RLIMIT_AS, &rl) == -1)
                    printf("WARNING! Could not set resource limit: Virtual memory.\n");
            } }

         const char *error;
         Solver S;
         S.max_decision_var = opt_restrict_decisions;
#ifdef DEBUG_SOLVER
         S.dbg_solver = new Solver();
#endif
         gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
             if (in == NULL)
                 printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);

             if (S.verbosity > 0){
                 printf("============================[ Problem Statistics ]=============================\n");
                 printf("|                                                                             |\n"); }

             parse_GRAPH(in, S);
             gzclose(in);

             if(strlen(graphstr)){
				 gzFile gin =gzopen(graphstr, "rb");
				parse_GRAPH(gin,S);
				gzclose(gin);
             }

         // Change to signal-handlers that will only notify the solver and allow it to terminate
           // voluntarily:
           signal(SIGINT, SIGINT_interrupt);
           signal(SIGXCPU,SIGINT_interrupt);
     //      printf("Solving circuit with %d gates, %d latches, %d inputs, %d outputs\n", aiger->num_ands, aiger->num_latches, aiger->num_inputs, aiger->num_outputs);

           const char * priority_file = opt_priority;
           if(strlen(priority_file)>0){
        	   FILE * f = fopen(priority_file,"r");
        	   if(f){
        		   char * line = NULL;
				  int v=0;
				  int p=0;
				  int total_read = 0;
				  while (fscanf(f," %d %d ", &v,&p) ==2) {
					  if(v<1 || v> S.nVars() || p<0){
						  fprintf(stderr,"Bad priority line: %d %d", v, p);
						  exit(1);
					  }
					  v--;
					  total_read++;
					  S.setDecisionPriority(v,p);
					 }
				  if(total_read==0){
					  fprintf(stderr,"Warning: read no priorities from priority file!");
				  }
        		   fclose(f);
        	   }else{
        		   fprintf(stderr,"Failed to read priority file!\n");
        	   }
           }


         //really simple, unsophisticated incremental BMC:

           lbool ret=S.solve()?l_True:l_False;


        if(ret==l_True){
        	printf("SAT\n");

        	if(S.theories.size()){
				Theory * t = S.theories[0];
				GraphTheorySolver *g = (GraphTheorySolver*)t;
				int width = sqrt(g->nNodes());
				int v = 0;
				//for (int i = 0;i<w;i++){
				//	for(int j = 0;j<w;j++){
				int lasty= 0;
				for(int n = 0;n<g->nNodes();n++){
					int x = n%width;
					int y = n/width;
					if(y > lasty)
						printf("\n");
						if (isatty(fileno(stdout))){

							if(S.model[n]==l_True)
								printf("\033[1;42m\033[1;37m 1\033[0m");
							else
								printf("\033[1;44m\033[1;37m 0\033[0m");
						}else{

							if(S.model[n]==l_True)
								printf(" 1");
							else
								printf(" 0");
						}

					lasty=y;
				}
				printf("\n\n");
				if(opt_check_solution){
					if(!g->check_solved()){
						fprintf(stderr,"Error! Solution doesn't satisfy graph properties!\n");
						exit(1);
					}
				}

				if(opt_print_reach){
				 v = 0;
				//for (int i = 0;i<w;i++){
				//	for(int j = 0;j<w;j++){
				 lasty= 0;
				for(int n = 0;n<g->nNodes();n++){
					int x = n%width;
					int y = n/width;
					if(y > lasty)
						printf("\n");
						if (isatty(fileno(stdout))){

							if(S.model[n]==l_True)
								printf("\033[1;42m\033[1;37m%3d\033[0m",n);
							else
								printf("\033[1;44m\033[1;37m%3d\033[0m",n);
						}else{

							if(S.model[n]==l_True)
								printf(" 1");
							else
								printf(" 0");
						}

					lasty=y;
				}
				printf("\n");

				printf("\n");
				for(int t = 0;t<S.theories.size();t++){
					printf("Theory %d\n", t);
					GraphTheorySolver *g = (GraphTheorySolver*)S.theories[t];

					for(int r = 0;r<g->reach_detectors.size();r++){

						int width = sqrt(g->nNodes());
						int lasty= 0;
						int extra =  g->nNodes() % width ? (width- g->nNodes() % width ):0;
						for(int n = 0;n<g->nNodes();n++){
							int x = n%width;

							int y = (n + extra )/width;
							if(y > lasty)
								printf("\n");

							int v =var( g->reach_detectors[r]->reach_lits[n]);
							if (isatty(fileno(stdout))){
								if(S.model[v]==l_True)
									printf("\033[1;42m\033[1;37m%4d\033[0m", v+1);
								else
									printf("\033[1;44m\033[1;37m%4d\033[0m",v+1);
							}else{

								if(S.model[v]==l_True)
									printf(" 1");
								else
									printf(" 0");
							}

							lasty=y;
						}
						printf("\n");
					}



					//g->drawFull();

					assert(g->dbg_solved());
				}

				}
/*        		for(int r = 0;r<g->reach_detectors.size();r++){

					int width = sqrt(g->nNodes());
					int lasty= 0;
					int extra =  g->nNodes() % width ? (width- g->nNodes() % width ):0;
					for(int n = 0;n<g->nNodes();n++){
						int x = n%width;

						int y = (n + extra )/width;

						int v =var( g->reach_detectors[r]->reach_lits[n]);
						if(v==306){
							int a =1;
						}
						if(S.model[v]==l_True){
							assert(S.value(v)==l_True);
							int node = g->reach_detectors[r]->getNode(v);
							g->reach_detectors[r]->positive_reach_detector->dbg_path(node);
							int  b=1;

						}


						lasty=y;
					}
					printf("\n");
				}*/
        	}

        }else if(ret==l_False){
        	printf("UNSAT\n");
        }else{
        	printf("UNKNOWN\n");
        }
        printStats(S);
        for(int i = 0;i<S.theories.size();i++){
        	Theory * t = S.theories[i];
        	GraphTheorySolver *g = (GraphTheorySolver*)t;
        	g->printStats();
        }
        fflush(stdout);
#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}
