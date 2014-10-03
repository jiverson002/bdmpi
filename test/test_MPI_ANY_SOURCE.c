#include <GKlib.h>
#include <math.h>
#include <stdio.h>
#include <bdmpi.h>

/************************************************************
 * THIS CODE WAS RUNNED USING THE FOLLOWING COMMAND LINE
 bdmprun -ns 4 build/Linux-x86_64/test/test_MPI_ANY_SOURCE
 ***************************************************************/

float fct(float x)
{
      return cos(x);
}

/* Prototype */
void other_work(char* header);
float integral(float a, int i, float h, int n);

int main(int argc, char* argv[])
{
      int n, p, myid, tag, master, proc;
      float h, integral_sum, a, b, pi, my_int;
      MPI_Request req;
      MPI_Status status;

/* Starts MPI processes ... */

      MPI_Init(&argc,&argv);                 /* starts MPI */
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);  /* get current process id */
      MPI_Comm_size(MPI_COMM_WORLD, &p);     /* get number of processes */

      master = 0;
      pi = acos(-1.0);  /* = 3.14159... */
      a = 0.;           /* lower limit of integration */
      b = pi*1./2.;     /* upper limit of integration */
      n = 500;          /* number of increment within each process */
      tag = 123;        /* set the tag to identify this particular job */
      h = (b-a)/n/p;    /* length of increment */

      my_int = integral(a,myid,h,n);  /* 0<=myid<=p-1 */

      printf("Process %d has the partial result of %f\n", myid, my_int);

      if(myid == master) {
        integral_sum = my_int;
        for (proc=1;proc<p;proc++) {
          MPI_Recv(&my_int, 1, MPI_FLOAT, MPI_ANY_SOURCE, tag,
                   MPI_COMM_WORLD, &status);
          printf("[%d] Message %f received %zu from %d\n", myid, my_int,
            status.count, status.MPI_SOURCE);
          /* with MPI_ANY_SOURCE, more efficient and less prone to deadlock */
          integral_sum += my_int;
        }
        printf("The Integral =%f\n",integral_sum);
      }
      else {
        MPI_Isend(&my_int, 1, MPI_FLOAT, master, tag,
                      MPI_COMM_WORLD, &req);      /* send my_int to master */
	      printf("[%d] Sending...\n", myid);
        MPI_Wait(&req, &status);
      }
      MPI_Finalize();                       /* let MPI finish up ... */
}

float integral(float a, int i, float h, int n)
{
      int j;
      float h2, aij, integ;

      integ = 0.0;                 /* initialize integral */
      h2 = h/2.;
      for (j=0;j<n;j++) {             /* sum over all "j" integrals */
        aij = a + (i*n + j)*h;        /* lower limit of "j" integral */
        integ += fct(aij+h2)*h;
      }
      return integ;
}
