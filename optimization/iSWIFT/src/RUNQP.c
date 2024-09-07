#include "Prime.h"
#include "Matrices.h"

int main(){
	double tset = 0.0, tsol = 0.0, t_ldl = 0.0;
	
		QP* myQP;
		myQP = QP_SETUP(n, m, p, Pjc, Pir, Ppr, Ajc, Air, Apr, Gjc, Gir, Gpr, c, h, b, sigma_d, P);

		idxint ExitCode = QP_SOLVE(myQP);
		if (myQP != NULL)
			printf("Setup Time     : %f ms\n", myQP->stats->tsetup*1000.0);
		if (ExitCode == QP_OPTIMAL){
			printf("Solve Time     : %f ms\n", (myQP->stats->tsolve + myQP->stats->tsetup)*1000.0);
			printf("KKT_Solve Time : %f ms\n", myQP->stats->kkt_time*1000.0);
			printf("LDL Time       : %f ms\n", myQP->stats->ldl_numeric*1000.0);
			printf("Diff	       : %f ms\n", (myQP->stats->kkt_time - myQP->stats->ldl_numeric)*1000.0);
			printf("Iterations     : %d\n", myQP->stats->IterationCount);
			printf("Optimal Solution Found\n");
		}
		if (ExitCode == QP_MAXIT){
			printf("Solve Time     : %f ms\n", myQP->stats->tsolve*1000.0);
			printf("KKT_Solve Time : %f ms\n", myQP->stats->kkt_time*1000.0);
			printf("LDL Time       : %f ms\n", myQP->stats->ldl_numeric*1000.0);
			printf("Diff	       : %f ms\n", (myQP->stats->kkt_time - myQP->stats->ldl_numeric)*1000.0);
			printf("Iterations     : %d\n", myQP->stats->IterationCount);
			printf("Maximum Iterations reached\n");
		}

		if (ExitCode == QP_FATAL){
			printf("Unknown Error Detected\n");
		}

		if (ExitCode == QP_KKTFAIL){
			printf("LDL Factorization fail\n");
		}



		QP_CLEANUP(myQP);

	

}
